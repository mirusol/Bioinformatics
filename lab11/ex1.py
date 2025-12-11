# moodle exercise

from __future__ import annotations
import argparse
from dataclasses import dataclass
from typing import Iterable, List, Tuple
@dataclass
class AlignmentResult:
	aligned_seq1: str
	aligned_seq2: str
	matches: int
	length: int
	similarity: float
	score_matrix: List[List[int]]
	path: List[Tuple[int, int]]


def needleman_wunsch(
	seq1: str,
	seq2: str,
	match_score: int = 1,
	mismatch_score: int = -1,
	gap_penalty: int = 0,
) -> AlignmentResult:

	n, m = len(seq1), len(seq2)
	scores: List[List[int]] = [[0] * (m + 1) for _ in range(n + 1)]
	pointers: List[List[str]] = [[""] * (m + 1) for _ in range(n + 1)]

	for i in range(1, n + 1):
		scores[i][0] = scores[i - 1][0] + gap_penalty
		pointers[i][0] = "up"
	for j in range(1, m + 1):
		scores[0][j] = scores[0][j - 1] + gap_penalty
		pointers[0][j] = "left"

	for i in range(1, n + 1):
		for j in range(1, m + 1):
			diag = scores[i - 1][j - 1] + (
				match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score
			)
			up = scores[i - 1][j] + gap_penalty
			left = scores[i][j - 1] + gap_penalty
			best = max(diag, up, left)

			if best == diag:
				pointers[i][j] = "diag"
			elif best == up:
				pointers[i][j] = "up"
			else:
				pointers[i][j] = "left"

			scores[i][j] = best

	aligned1: List[str] = []
	aligned2: List[str] = []
	path: List[Tuple[int, int]] = []
	i, j = n, m
	while i > 0 or j > 0:
		path.append((i, j))
		move = pointers[i][j]
		if move == "diag":
			aligned1.append(seq1[i - 1])
			aligned2.append(seq2[j - 1])
			i -= 1
			j -= 1
		elif move == "up":
			aligned1.append(seq1[i - 1])
			aligned2.append("-")
			i -= 1
		else:  # left
			aligned1.append("-")
			aligned2.append(seq2[j - 1])
			j -= 1

	path.append((0, 0))
	path.reverse()
	aligned_seq1 = "".join(reversed(aligned1))
	aligned_seq2 = "".join(reversed(aligned2))

	matches = sum(a == b for a, b in zip(aligned_seq1, aligned_seq2))
	length = len(aligned_seq1)
	similarity = 100.0 * matches / length if length else 0.0

	return AlignmentResult(
		aligned_seq1=aligned_seq1,
		aligned_seq2=aligned_seq2,
		matches=matches,
		length=length,
		similarity=similarity,
		score_matrix=scores,
		path=path,
	)


def _match_bar(row1: str, row2: str) -> str:
	return "".join("|" if a == b else " " for a, b in zip(row1, row2))


def _plot_alignment(result: AlignmentResult, seq1: str, seq2: str) -> None:
	try:
		import matplotlib.pyplot as plt
		from matplotlib.colors import ListedColormap
		import numpy as np
	except ImportError:
		print("matplotlib not available; skipping plots.")
		return

	scores = np.array(result.score_matrix)
	path = set(result.path)
	path_matrix = np.zeros_like(scores)
	for (i, j) in path:
		path_matrix[i, j] = 1

	fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 5))

	im0 = ax0.imshow(scores, cmap="magma", origin="upper")
	ax0.set_title("Score matrix")
	ax0.set_xticks(range(len(seq2) + 1))
	ax0.set_yticks(range(len(seq1) + 1))
	ax0.set_xticklabels([" "] + list(seq2))
	ax0.set_yticklabels([" "] + list(seq1))
	ax0.tick_params(axis="x", labelrotation=90)
	fig.colorbar(im0, ax=ax0, fraction=0.046, pad=0.04)
	
	cmap = ListedColormap([[1.0, 0.98, 0.85], [0.8, 0.0, 0.0]])
	im1 = ax1.imshow(path_matrix, cmap=cmap, origin="upper")
	ax1.set_title("Traceback path")
	ax1.set_xticks(range(len(seq2) + 1))
	ax1.set_yticks(range(len(seq1) + 1))
	ax1.set_xticklabels(range(len(seq2) + 1))
	ax1.set_yticklabels(range(len(seq1) + 1))
	ax1.grid(which="both", color="black", linewidth=0.5)
	ax1.set_xlim(-0.5, len(seq2) + 0.5)
	ax1.set_ylim(len(seq1) + 0.5, -0.5)

	fig.tight_layout()
	plt.show()


def format_alignment(result: AlignmentResult) -> str:
	bar = _match_bar(result.aligned_seq1, result.aligned_seq2)
	lines = [
		result.aligned_seq1,
		bar,
		result.aligned_seq2,
		f"Matches = {result.matches}",
		f"Length = {result.length}",
		f"Similarity = {result.similarity:.0f} %",
	]
	return "\n".join(lines)


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
	parser = argparse.ArgumentParser(description="Needleman-Wunsch alignment")
	parser.add_argument("--seq1", default="ACCGTGAAGCCAATAC", help="First DNA sequence")
	parser.add_argument("--seq2", default="AGCGTGCAGCCAATAC", help="Second DNA sequence")
	parser.add_argument("--match", type=int, default=1, help="Score for a match")
	parser.add_argument("--mismatch", type=int, default=-1, help="Score for a mismatch")
	parser.add_argument("--gap", type=int, default=0, help="Gap penalty")
	parser.add_argument("--no-plot", action="store_true", help="Skip matplotlib plots")
	return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> None:
	args = parse_args(argv)
	result = needleman_wunsch(
		args.seq1,
		args.seq2,
		match_score=args.match,
		mismatch_score=args.mismatch,
		gap_penalty=args.gap,
	)

	print("Alignment:\n")
	print(format_alignment(result))

	if not args.no_plot:
		_plot_alignment(result, args.seq1, args.seq2)


if __name__ == "__main__":
	main()

