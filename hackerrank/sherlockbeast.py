n = int(raw_input().strip())
for i in range(n):
    seq = list(raw_input().strip())
    out = [seq[0]]
    for j in range(1, len(seq)):
        if seq[j] != out[-1]:
            out.append(seq[j])
    print len(seq)-len(out)