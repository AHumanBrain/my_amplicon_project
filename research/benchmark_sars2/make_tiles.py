"""Overlapping 400 bp windows (step 300) over SARS-CoV-2 NC_045512.2: python make_tiles.py > sars2_tiles.tsv"""
L = 29903
print("name	reference	start	end	amplicon_min	amplicon_max")
for i in range(100):
    s = 1 + i * 300
    e = min(s + 399, L)
    if e - s < 300:
        break
    print(f"W{i + 1:03d}	NC_045512.2	{s}	{e}	330	400")
