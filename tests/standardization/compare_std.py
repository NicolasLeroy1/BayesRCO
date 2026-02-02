import sys

def parse_std_output(filename):
    sections = {}
    current_section = None
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.endswith(':'):
                current_section = line[:-1].strip()
                sections[current_section] = []
            elif current_section:
                try:
                    sections[current_section].append(float(line))
                except ValueError:
                    pass
    return sections

def compare(f1, f2):
    s_f = parse_std_output(f1)
    s_c = parse_std_output(f2)

    if set(s_f.keys()) != set(s_c.keys()):
        print("ERROR: Section mismatch!")
        print(f"Fortran sections: {list(s_f.keys())}")
        print(f"C sections: {list(s_c.keys())}")
        return False

    all_match = True
    for section in s_f:
        v_f = s_f[section]
        v_c = s_c[section]

        if len(v_f) != len(v_c):
            print(f"ERROR: Length mismatch in section {section}!")
            all_match = False
            continue

        for i, (val1, val2) in enumerate(zip(v_f, v_c)):
            if abs(val1 - val2) > 1e-12:
                print(f"ERROR: Mismatch in {section} at index {i}: Fortran={val1:.16f}, C={val2:.16f}")
                all_match = False

    if all_match:
        print("SUCCESS: Standardization is equivalent between Fortran and C!")
    return all_match

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python compare_std.py std_f.txt std_c.txt")
        sys.exit(1)
    if compare(sys.argv[1], sys.argv[2]):
        sys.exit(0)
    else:
        sys.exit(1)
