import struct
import random

def generate_plink_data(prefix, n_ind=10, n_loci=20):
    # .fam file: FID IID PID MID SEX PHENOTYPE
    with open(f"{prefix}.fam", "w") as f:
        for i in range(1, n_ind + 1):
            pheno_choice = random.choice([1.5, 2.5, -9.0, None])
            pheno_str = "NA" if pheno_choice is None else f"{pheno_choice:.2f}"
            f.write(f"FAM001 IND{i:03d} 0 0 1 {pheno_str}\n")

    # .bim file: CHR ID POS GDIST A1 A2
    with open(f"{prefix}.bim", "w") as f:
        for j in range(1, n_loci + 1):
            f.write(f"1 SNP{j:03d} 0 {j*100} A G\n")

    # .bed file: Magic (3 bytes) + packed genotypes
    with open(f"{prefix}.bed", "wb") as f:
        f.write(bytes([0x6c, 0x1b, 0x01])) # Magic number snp-major
        
        mapping = {0: 0b00, 1: 0b10, 2: 0b11, 3: 0b01}
        
        for j in range(n_loci):
            genotypes = [random.randint(0, 3) for _ in range(n_ind)]
            codes = [mapping[g] for g in genotypes]
            
            n_bytes = (n_ind + 3) // 4
            for b_idx in range(n_bytes):
                byte_val = 0
                for i in range(4):
                    idx = b_idx * 4 + i
                    if idx < n_ind:
                        byte_val |= (codes[idx] << (2 * i))
                f.write(struct.pack("B", byte_val))

if __name__ == "__main__":
    random.seed(42)
    generate_plink_data("test_data")
    print("Generated test_data.bed, test_data.bim, test_data.fam")
