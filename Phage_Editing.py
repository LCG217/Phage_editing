import primer3
from Bio import SeqIO
import pandas as pd

# file for main sequence and any insert sequence
import primer3
from Bio import SeqIO

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

print("\n✅ Settings summary:")
print(f"Main sequence file: {file_path}")
print(f"Insert sequence: {insert_file if insert_file else 'Not used'}")
print(f"Exclusion zones: {exclusion_zones if exclusion_zones else 'None'}")


# define exclusion zones
exclusion_zones = [(34263, 36284)]
# overhang bases
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

def overlaps_exclusion(pos, exclusion_zones):
    """Check if primer position overlaps exclusion zone."""
    for start, end in exclusion_zones:
        if start <= pos <= end:
            return True
    return False

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

    # Handle exclusion zones and insert sequence logic
    if use_exclusions == "y" and (overlaps_exclusion(left_loc, exclusion_zones) or overlaps_exclusion(right_loc, exclusion_zones)):
        if use_insert == "y" and insert_file:
            insert_seq = get_insert_seq(insert_file)
            oh_start, oh_end = make_overhangs(insert_seq, flank=FLANK)
            forward = oh_start + forward
            reverse = oh_end + reverse
            print(f"Primer {primer_index}: Overlaps exclusion zone → overhangs added.")
        else:
            print(f"Primer {primer_index}: Overlaps exclusion zone but no insert sequence provided.")
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

output_str = ''
for key,value in all_primer_pairs.items():
    output_str += f">{key}_fwd\n{value['forward']}\n>{key}_rev\n{value['reverse']}\n"

print(output_str)