import primer3
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import random

# Prompt the user for file paths instead of hardcoding
file_path = input("Enter the path to your main sequence file (e.g. genbank or FASTA): ").strip()
use_insert = input("Do you want to include an insert sequence? (y/n): ").strip().lower()
insert_file = None
if use_insert == "y":
    insert_file = input("Enter the path to your insert sequence file in FASTA format: ")

# exclusion zone selection
use_exclusions = input("Do you want to include exclusion zones? (y/n): ").strip().lower()
exclusion_zones = []
if use_exclusions == "y":
    print("Enter exclusion zones as start:end (e.g., 34263:36284). Enter 'done' when finished.")
    while True:
        zone = input("Exclusion zone: ").strip()
        if zone.lower() == "done":
            break
        try:
            start, end = map(int, zone.split(":"))
            exclusion_zones.append((start, end))
        except ValueError:
            print("Invalid format. Use start:end (e.g., 100:500).")

print("Settings summary:")
print(f"Main sequence file: {file_path}")
print(f"Insert sequence: {insert_file if insert_file else 'Not used'}")
print(f"Exclusion zones: {exclusion_zones if exclusion_zones else 'None'}")

# number of overhang bases
FLANK = 20

# Primer design parameters
Opt_Size = 24
Min_Size = 18
Max_Size = 35
Opt_TM = 60
Min_TM = 55
Max_TM = 65
Min_GC = 40
Max_GC = 60
size_range_min = 4500
size_range_max = 6000
force_end = -1000000
force_primers = 1

all_primer_pairs = {}
idt_primers = []
# index for primer naming _1, _2 etc
primer_index = 1

# Dictionary to store random overhangs for each exclusion zone
exclusion_overhangs = {}

# read sequences, define genome length and get sequences in the right format
record = SeqIO.read(file_path, 'gb')
Seq_Template = str(record.seq)
genomelength = len(Seq_Template)
print(f"Genome length: {genomelength}")
nuc_region_of_interest = record.seq


# define functions
def get_insert_seq(insert_file):
    insert_rec = SeqIO.read(insert_file, 'fasta')
    return insert_rec.seq


def make_overhangs(insert_seq, flank=20):
    """Return forward/reverse overhangs from insert ends."""
    start_overhang = insert_seq[:flank].reverse_complement()
    end_overhang = insert_seq[-flank:]
    return str(start_overhang), str(end_overhang)


def generate_random_overhang(length=20):
    """Generate a random DNA sequence of given length."""
    return ''.join(random.choice('ATGC') for _ in range(length))


def make_complementary_overhangs_for_gibson(Seq_Template, left_loc, right_loc, flank=20):
    # Forward primer overhang = upstream region before left_loc
    start_oh_start = max(0, left_loc - flank)
    fwd_overhang = Seq_Template[start_oh_start:left_loc]

    # Reverse primer overhang = downstream region after right_loc, reverse-complemented
    end_oh_end = min(len(Seq_Template), right_loc + flank)
    rev_overhang_raw = Seq_Template[right_loc:end_oh_end]
    rev_overhang = str(Seq(rev_overhang_raw).reverse_complement())

    return fwd_overhang, rev_overhang


def overlaps_exclusion(pos, exclusion_zones):
    """Check if primer position overlaps exclusion zone."""
    for start, end in exclusion_zones:
        if start <= pos <= end:
            return True
    return False


def get_spanning_exclusion_zone(left_loc, right_loc, exclusion_zones):
    """
    Check if this amplicon spans across an exclusion zone.
    Returns the exclusion zone if:
    - Left primer is before the exclusion start, AND
    - Right primer is at or after the exclusion start (but could be before or after the end)
    This identifies the LAST primer pair before the deletion.
    """
    for start, end in exclusion_zones:
        if left_loc < start and right_loc >= start:
            return (start, end)
    return None


def starts_after_exclusion_zone(left_loc, right_loc, exclusion_zones):
    """
    Check if this amplicon starts right after an exclusion zone.
    Returns the exclusion zone if:
    - Left primer is at or just after the exclusion end
    This identifies the FIRST primer pair after the deletion.
    """
    for start, end in exclusion_zones:
        # Left primer should be at or very close to the end of the exclusion
        if left_loc >= end and left_loc - end <= 100:  # Within 100bp of exclusion end
            return (start, end)
    return None


def primer_location_start(location):
    for start, end in exclusion_zones:
        if (start - 65) <= location < end:
            return end
    return max(0, location - 65)


def primer_location_end(size_range_min, Seq_Start, force_end):
    for start, end in exclusion_zones:
        if Seq_Start < start:
            if Seq_Start + size_range_min > start:
                return start - 1
        elif end < Seq_Start < start:
            if Seq_Start + size_range_min > start:
                return start - 1
    return force_end


def size_range(size_range_min, Seq_Start):
    for start, end in exclusion_zones:
        if Seq_Start < start:
            if Seq_Start + size_range_min > start:
                return start - Seq_Start
        elif end < Seq_Start < start:
            if Seq_Start + size_range_min > start:
                return start - Seq_Start
    return size_range_min


# function to create excel file of primers for easy importing/ordering
def export_primers_to_excel(primer_dict, filename="primers.xlsx"):
    rows = []

    for name, pair in primer_dict.items():
        rows.append({"Primer Name": f"{name}_F", "Primer Sequence": pair['forward']})
        rows.append({"Primer Name": f"{name}_R", "Primer Sequence": pair['reverse']})

    df = pd.DataFrame(rows)
    df.to_excel(filename, index=False)

    print(f"Primers saved to Excel file: {filename}")


# loop for main code
current_location = 0

while current_location < genomelength:
    # Set starting point for this round
    Seq_Start = primer_location_start(current_location) if use_exclusions == "y" else max(0, current_location - 65)
    force_end = primer_location_end(size_range_min, Seq_Start, force_end) if use_exclusions == "y" else force_end
    size_range_min = size_range(size_range_min, Seq_Start) if use_exclusions == "y" else size_range_min
    Size_Range = (size_range_min, size_range_max)

    # Update SEQUENCE_INCLUDED_REGION for the current window
    include_length = genomelength - Seq_Start
    if include_length < size_range_min:
        print(f"Remaining sequence too short (<{size_range_min}bp) at position {Seq_Start}. Stopping.")
        break
    Seq_Inc_Region = (Seq_Start, include_length)

    # Primer3 input
    Seq_ID = f'T7_{primer_index}'
    seq_args = {
        'SEQUENCE_ID': Seq_ID,
        'SEQUENCE_TEMPLATE': str(nuc_region_of_interest),
        'SEQUENCE_INCLUDED_REGION': Seq_Inc_Region,
        'SEQUENCE_FORCE_LEFT_START': Seq_Start,
        'SEQUENCE_FORCE_RIGHT_END': force_end,
    }

    global_args = {
        'PRIMER_OPT_SIZE': Opt_Size,
        'PRIMER_MIN_SIZE': Min_Size,
        'PRIMER_MAX_SIZE': Max_Size,
        'PRIMER_OPT_TM': Opt_TM,
        'PRIMER_MIN_TM': Min_TM,
        'PRIMER_MAX_TM': Max_TM,
        'PRIMER_MIN_GC': Min_GC,
        'PRIMER_MAX_GC': Max_GC,
        'PRIMER_PRODUCT_SIZE_RANGE': Size_Range,
        'PRIMER_PICK_ANYWAY': force_primers,
    }

    # Run primer3
    try:
        primers = primer3.bindings.design_primers(seq_args, global_args)
        forward = primers['PRIMER_LEFT_0_SEQUENCE']
        reverse = primers['PRIMER_RIGHT_0_SEQUENCE']
        left_loc = primers['PRIMER_LEFT_0'][0]
        right_loc = primers['PRIMER_RIGHT_0'][0]
    except KeyError:
        print(f"Primer {primer_index}: Primer design failed or no primer found. Stopping.")
        break

    # Handle exclusion zones and insert/overhang logic
    if use_exclusions == "y":
        # Check if this amplicon spans INTO an exclusion zone (last primer before deletion)
        spanning_zone = get_spanning_exclusion_zone(left_loc, right_loc, exclusion_zones)
        # Check if this amplicon starts right AFTER an exclusion zone (first primer after deletion)
        after_zone = starts_after_exclusion_zone(left_loc, right_loc, exclusion_zones)

        if spanning_zone and use_insert == "y" and insert_file:
            # Use insert sequence overhangs
            insert_seq = get_insert_seq(insert_file)
            oh_start, oh_end = make_overhangs(insert_seq, flank=FLANK)
            reverse = oh_end + reverse
            exclusion_overhangs[spanning_zone] = oh_end  # Store for next primer
            print(f"Primer {primer_index}: Spans into exclusion zone → insert overhang added to reverse primer.")

        elif spanning_zone and use_insert == "n":
            # Generate random overhang for reverse primer before exclusion
            random_overhang = generate_random_overhang(FLANK)
            reverse = random_overhang + reverse
            exclusion_overhangs[spanning_zone] = random_overhang  # Store for next primer
            print(
                f"Primer {primer_index}: Spans into exclusion zone → random overhang added to reverse primer: {random_overhang}")

        elif after_zone:
            # Add reverse complement of stored overhang to forward primer
            if after_zone in exclusion_overhangs:
                stored_overhang = exclusion_overhangs[after_zone]
                rc_overhang = str(Seq(stored_overhang).reverse_complement())
                forward = rc_overhang + forward
                print(
                    f"Primer {primer_index}: After exclusion zone → reverse complement overhang added to forward primer: {rc_overhang}")
            else:
                print(f"Primer {primer_index}: After exclusion zone but no stored overhang found.")
        else:
            print(f"Primer {primer_index}: {forward} / {reverse} at {right_loc}")
    else:
        print(f"Primer {primer_index}: {forward} / {reverse} at {right_loc}")

    # Store results
    name = f"primer_pair_{primer_index}"
    all_primer_pairs[name] = {'forward': forward, 'reverse': reverse}

    idt_primers.append([f'Forward Primer {primer_index}', forward, '25nm', 'STD'])
    idt_primers.append([f'Reverse Primer {primer_index}', reverse, '25nm', 'STD'])

    # Prepare for next loop
    primer_index += 1
    current_location = right_loc - 65

export_primers_to_excel(all_primer_pairs, "primer_output.xlsx")


