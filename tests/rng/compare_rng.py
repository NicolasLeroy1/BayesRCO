import sys

def parse_output(filename):
    sections = {}
    current_section = None
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.endswith(':'):
                current_section = line[:-1]
                sections[current_section] = []
            elif current_section:
                try:
                    sections[current_section].append(float(line))
                except ValueError:
                    pass
    return sections

def compare(f1, f2):
    sections_f = parse_output(f1)
    sections_c = parse_output(f2)

    if set(sections_f.keys()) != set(sections_c.keys()):
        print("ERROR: Section mismatch!")
        print(f"Fortran sections: {list(sections_f.keys())}")
        print(f"C sections: {list(sections_c.keys())}")
        return False

    all_match = True
    for section in sections_f:
        vals_f = sections_f[section]
        vals_c = sections_c[section]

        if len(vals_f) != len(vals_c):
            print(f"ERROR: Length mismatch in section {section}!")
            all_match = False
            continue

        for i, (v_f, v_c) in enumerate(zip(vals_f, vals_c)):
            if abs(v_f - v_c) > 1e-15:
                print(f"ERROR: Mismatch in {section} at index {i}: Fortran={v_f}, C={v_c}")
                all_match = False

    if all_match:
        print("SUCCESS: All RNG values match bit-precisely (or within machine epsilon)!")
    return all_match

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python compare_rng.py results_f.txt results_c.txt")
        sys.exit(1)
    
    if compare(sys.argv[1], sys.argv[2]):
        sys.exit(0)
    else:
        sys.exit(1)
