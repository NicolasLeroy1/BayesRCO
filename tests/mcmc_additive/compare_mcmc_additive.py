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
                current_section = line[:-1].strip()
                sections[current_section] = []
            elif current_section:
                try:
                    sections[current_section].append(float(line))
                except ValueError:
                    try:
                        sections[current_section].append(int(line))
                    except ValueError:
                        pass
    return sections

def compare(f1, f2, tol=1e-10):
    s_f = parse_output(f1)
    s_c = parse_output(f2)

    if set(s_f.keys()) != set(s_c.keys()):
        print("ERROR: Section mismatch!")
        print(f"Fortran sections: {list(s_f.keys())}")
        print(f"C sections: {list(s_c.keys())}")
        return False

    all_match = True
    for section in sorted(s_f.keys()):
        v_f = s_f[section]
        v_c = s_c[section]

        if len(v_f) != len(v_c):
            print(f"ERROR: Length mismatch in section {section}! (F={len(v_f)}, C={len(v_c)})")
            all_match = False
            continue

        section_ok = True
        for i, (val1, val2) in enumerate(zip(v_f, v_c)):
            diff = abs(val1 - val2)
            if diff > tol:
                print(f"MISMATCH in {section}[{i}]: Fortran={val1}, C={val2}, diff={diff:.2e}")
                all_match = False
                section_ok = False

        if section_ok:
            print(f"OK: {section} ({len(v_f)} values match)")

    if all_match:
        print("\nSUCCESS: MCMC additive kernel is equivalent between Fortran and C!")
    else:
        print("\nFAILURE: Differences detected!")
    return all_match

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python compare_mcmc_additive.py results_f.txt results_c.txt")
        sys.exit(1)
    if compare(sys.argv[1], sys.argv[2]):
        sys.exit(0)
    else:
        sys.exit(1)
