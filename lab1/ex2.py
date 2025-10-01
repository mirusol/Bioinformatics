s = "ACGGGCATATGCGC".replace(" ", "")

alphabet = []
for c in s:
    if c not in alphabet:
        alphabet.append(c)

print("alphabet:", ''.join(alphabet))

for letter in alphabet:
    percent = s.count(letter) * 100 / len(s)
    print(f"{letter}: {percent:.2f}%")