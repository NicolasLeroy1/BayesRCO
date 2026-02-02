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

    # Known index offset adjustments:
    # Fortran uses 1-based indexing, C uses 0-based
    # a_after_init: Fortran non-zero values should be compared with C value + 1
    # (Or equivalently, Fortran value - 1 should equal C value for non-zero Fortran)
    # For uninitialized (both 0), they should match as-is
    index_offset_sections = {'a_after_init'}

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
            # For a_after_init: 
            # - Fortran stores 1-based index, C stores 0-based
            # - Fortran 0 means uninitialized (nannot > 1), C 0 also means uninitialized
            # - So: Fortran non-zero should equal C + 1, both zeros should match
            if section in index_offset_sections:
                if val1 == 0 and val2 == 0:
                    # Both uninitialized, match
                    continue
                elif val1 != 0:
                    # Fortran is 1-based, so Fortran - 1 should equal C
                    adjusted_val1 = val1 - 1
                    diff = abs(adjusted_val1 - val2)
                else:
                    # Fortran is 0 but C is non-zero - mismatch
                    diff = abs(val1 - val2)
                    adjusted_val1 = val1
            else:
                adjusted_val1 = val1
                diff = abs(adjusted_val1 - val2)
            
            if diff > tol:
                if section in index_offset_sections:
                    print(f"MISMATCH in {section}[{i}]: Fortran={val1} (adj={adjusted_val1}), C={val2}, diff={diff:.2e}")
                else:
                    print(f"MISMATCH in {section}[{i}]: Fortran={val1}, C={val2}, diff={diff:.2e}")
                all_match = False
                section_ok = False

        if section_ok:
            if section in index_offset_sections:
                print(f"OK: {section} ({len(v_f)} values match, adjusted for 1-based/0-based offset)")
            else:
                print(f"OK: {section} ({len(v_f)} values match)")

    if all_match:
        print("\nSUCCESS: MCMC mixture kernel is equivalent between Fortran and C!")
    else:
        print("\nFAILURE: Differences detected!")
    return all_match

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python compare_mcmc_mixture.py results_f.txt results_c.txt")
        sys.exit(1)
    if compare(sys.argv[1], sys.argv[2]):
        sys.exit(0)
    else:
        sys.exit(1)
