# 
from __future__ import annotations
import argparse
import os
from dataclasses import dataclass
from typing import Iterable, List, Tuple


@dataclass
class LocalAlignment:
	aligned_a: str
	aligned_b: str
	score: int
	start_a: int
	end_a: int
	start_b: int
	end_b: int
	identity: float
	penalty_similarity: float
	norm_score: float
	coverage: float


def read_fasta(path: str) -> str:
	with open(path, "r", encoding="utf-8") as f:
		lines = [line.strip() for line in f if line.strip() and not line.startswith(">")]
	return "".join(lines).upper()


def smith_waterman(
	a: str,
	b: str,
	match: int = 2,
	mismatch: int = -1,
	gap: int = -2,
) -> LocalAlignment:
	na, nb = len(a), len(b)
	scores = [[0] * (nb + 1) for _ in range(na + 1)]
	pointers: List[List[int]] = [[0] * (nb + 1) for _ in range(na + 1)]

	best_score = 0
	best_pos = (0, 0)

	for i in range(1, na + 1):
		for j in range(1, nb + 1):
			diag = scores[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch)
			up = scores[i - 1][j] + gap
			left = scores[i][j - 1] + gap
			best = max(0, diag, up, left)
			scores[i][j] = best
			if best == 0:
				pointers[i][j] = 0
			elif best == diag:
				pointers[i][j] = 1
			elif best == up:
				pointers[i][j] = 2
			else:
				pointers[i][j] = 3
			if best > best_score:
				best_score = best
				best_pos = (i, j)

	i, j = best_pos
	aligned_a: List[str] = []
	aligned_b: List[str] = []
	end_a, end_b = i, j
	while i > 0 and j > 0 and scores[i][j] > 0:
		move = pointers[i][j]
		if move == 1:
			aligned_a.append(a[i - 1])
			aligned_b.append(b[j - 1])
			i -= 1
			j -= 1
		elif move == 2:
			aligned_a.append(a[i - 1])
			aligned_b.append("-")
			i -= 1
		else:
			aligned_a.append("-")
			aligned_b.append(b[j - 1])
			j -= 1

	start_a, start_b = i, j
	aligned_a.reverse()
	aligned_b.reverse()
	aligned_a_str = "".join(aligned_a)
	aligned_b_str = "".join(aligned_b)
	matches = sum(x == y and x != "-" and y != "-" for x, y in zip(aligned_a_str, aligned_b_str))
	mismatches = sum(x != y and x != "-" and y != "-" for x, y in zip(aligned_a_str, aligned_b_str))
	gaps = sum(x == "-" or y == "-" for x, y in zip(aligned_a_str, aligned_b_str))
	L = len(aligned_a_str)
	identity = 100.0 * matches / L if L else 0.0
	penalty_similarity = 100.0 * (matches - mismatches - gaps) / L if L else 0.0
	pos_score = matches * match + mismatches * mismatch + gaps * gap
	norm_score = 100.0 * pos_score / (L * max(match, 1)) if L else 0.0
	aligned_a_len = len([c for c in aligned_a_str if c != "-"])
	coverage = 100.0 * aligned_a_len / na if na else 0.0

	return LocalAlignment(
		aligned_a=aligned_a_str,
		aligned_b=aligned_b_str,
		score=best_score,
		start_a=start_a,
		end_a=end_a,
		start_b=start_b,
		end_b=end_b,
		identity=identity,
		penalty_similarity=penalty_similarity,
		norm_score=norm_score,
		coverage=coverage,
	)


def chunked_alignment(
	genome_a: str,
	genome_b: str,
	chunk_size: int,
	overlap: int,
	match: int,
	mismatch: int,
	gap: int,
) -> List[LocalAlignment]:
	if chunk_size <= 0 or overlap < 0 or overlap >= chunk_size:
		raise ValueError("chunk_size must be >0 and overlap in [0, chunk_size)")

	step = chunk_size - overlap
	results: List[LocalAlignment] = []
	max_pos = min(len(genome_a), len(genome_b))

	offset = 0
	while offset < max_pos:
		window_a = genome_a[offset : offset + chunk_size]
		window_b = genome_b[offset : offset + chunk_size]
		if not window_a or not window_b:
			break
		res = smith_waterman(window_a, window_b, match=match, mismatch=mismatch, gap=gap)
		res.start_a += offset
		res.end_a += offset
		res.start_b += offset
		res.end_b += offset
		results.append(res)
		offset += step
	return results


def plot_similarity(chunks: List[LocalAlignment], chunk_size: int, overlap: int, out_title: str) -> None:
	try:
		import matplotlib.pyplot as plt
	except ImportError:
		print("matplotlib not available; skipping plots.")
		return

	x = list(range(len(chunks)))
	y = [c.identity for c in chunks]
	starts = [c.start_a for c in chunks]
	ends = [c.end_a for c in chunks]

	fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
	ax0.plot(x, y, marker="o", linewidth=1.5)
	ax0.set_ylabel("Identity %")
	ax0.set_title(out_title)
	ax0.grid(True, linestyle=":", alpha=0.6)

	ax1.bar(x, [e - s for s, e in zip(starts, ends)], bottom=starts, width=0.8, color="#c44", alpha=0.7)
	ax1.set_ylabel("Influenza coord")
	ax1.set_xlabel("Chunk index (step = chunk_size - overlap)")
	ax1.grid(True, linestyle=":", alpha=0.4)
	plt.tight_layout()
	plt.show()


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
	parser = argparse.ArgumentParser(description="Chunked local alignment of influenza vs COVID-19 genomes")
	parser.add_argument("--influenza", default="influenza.fasta", help="Path to influenza genome FASTA")
	parser.add_argument("--covid", default="covid.fasta", help="Path to COVID-19 genome FASTA")
	parser.add_argument("--chunk-size", type=int, default=2000, help="Window size for stepwise alignment")
	parser.add_argument("--overlap", type=int, default=500, help="Overlap between consecutive windows")
	parser.add_argument("--match", type=int, default=2, help="Match score")
	parser.add_argument("--mismatch", type=int, default=-1, help="Mismatch score")
	parser.add_argument("--gap", type=int, default=-2, help="Gap penalty")
	parser.add_argument("--no-plot", action="store_true", help="Skip matplotlib plots")
	return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> None:
	args = parse_args(argv)

	if not os.path.exists(args.influenza) or not os.path.exists(args.covid):
		print("FASTA files not found. Place influenza.fasta and covid.fasta in this folder or use --influenza/--covid paths.")
		return

	influenza = read_fasta(args.influenza)
	covid = read_fasta(args.covid)

	print(f"Loaded influenza length: {len(influenza):,}")
	print(f"Loaded COVID-19 length:  {len(covid):,}")
	print(
		f"Running chunked local alignment with chunk={args.chunk_size}, overlap={args.overlap}, scores: match={args.match}, mismatch={args.mismatch}, gap={args.gap}"
	)

	chunks = chunked_alignment(
		influenza,
		covid,
		chunk_size=args.chunk_size,
		overlap=args.overlap,
		match=args.match,
		mismatch=args.mismatch,
		gap=args.gap,
	)

	if not chunks:
		print("No alignment chunks produced.")
		return

	print("\nTop 5 chunks by identity:")
	ranked = sorted(enumerate(chunks), key=lambda x: x[1].identity, reverse=True)[:5]
	for idx, res in ranked:
		print(
			f"Chunk {idx}: score={res.score}, identity={res.identity:.2f}%, penalty_sim={res.penalty_similarity:.2f}%, norm_score={res.norm_score:.2f}%, coverage={res.coverage:.2f}%, influenza[{res.start_a}:{res.end_a}], covid[{res.start_b}:{res.end_b}]"
		)

	if not args.no_plot:
		plot_similarity(chunks, args.chunk_size, args.overlap, "Influenza vs COVID-19 similarity by window")


if __name__ == "__main__":
	main()
