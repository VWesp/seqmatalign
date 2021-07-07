import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.ticker as ticker
import sys
import traceback

try:
    parser = argparse.ArgumentParser(description="Align two sequences with Needleman-Wunsch or Smith-Waterman and build the backtrace matrix")
    parser.add_argument("-s1", "--sequence1", help="Set the path to first sequence file (with or without header)", default="seq1.fasta")
    parser.add_argument("-s2", "--sequence2", help="Set the path to second sequence file (with or without header)", default="seq2.fasta")
    parser.add_argument("-alg", "--algorithm", help="Specify whether the Needleman-Wunsch (nw) or the Smith-Waterman (sw) algorithm should be used", choices=["nw", "sw"], default="nw")
    parser.add_argument("-f", "--function", help="Specify whether the score function should be minimized (distance, dis) or maximized (similarity, sim). If the Smith-Waterman algorithm is chosen, the function is always maximized", choices=["dis", "sim"], default="sim")
    parser.add_argument("-s", "--scores", help="Set the path to the score file", default="scores.txt")
    parser.add_argument("-g", "--gap", help="Specify the gap costs", default=-1, type=float)
    args = parser.parse_args()


    def highlight_cell(x,y, ax=None, **kwargs):
        rect = plt.Rectangle((x-.5, y-.5), 1,1, fill=False, **kwargs)
        ax = ax or plt.gca()
        ax.add_patch(rect)
        return rect


    seq1 = None
    with open(args.sequence1, "r") as seq_reader:
        content = seq_reader.readlines()
        if(content[0].startswith(">")):
            seq1 = "".join(content[1:]).strip().replace("\n", "")
        else:
            seq1 = "".join(content).strip().replace("\n", "")

    seq2 = None
    with open(args.sequence2, "r") as seq_reader:
        content = seq_reader.readlines()
        if(content[0].startswith(">")):
            seq2 = "".join(content[1:]).strip().replace("\n", "")
        else:
            seq2 = "".join(content).strip().replace("\n", "")

    scores = {}
    with open(args.scores, "r") as score_reader:
        content = score_reader.readlines()
        letters = content[0].strip().split(" ")
        first_index = 0
        for line in content[1:]:
            score = line.strip().split(" ")
            scores[letters[first_index]] = {}
            for second_index in range(len(score)):
                scores[letters[first_index]][letters[second_index]] = float(score[second_index])

            first_index += 1

    if(args.algorithm == "sw"):
        args.function = "sim"

    seq1 = "-" + seq1
    seq2 = "-" + seq2
    matrix = [[0 for i in range(len(seq2))] for i in range(len(seq1))]
    matrix_left = [[0 for i in range(len(seq2))] for i in range(len(seq1))]
    matrix_up = [[0 for i in range(len(seq2))] for i in range(len(seq1))]
    backtrace_matrix = [["" for i in range(len(seq2))] for i in range(len(seq1))]
    arrow_matrix = [[{} for i in range(len(seq2))] for i in range(len(seq1))]
    maximum_value = -(sys.maxsize - 1)
    maximum_indices = []
    figsize = max(len(seq1), len(seq2))/2
    fig, ax = plt.subplots(figsize=(figsize, figsize))
    ax.text(0, 0, str(matrix[0][0]), va='center', ha='center')
    for j in range(1, len(seq2)):
        if(args.algorithm == "nw"):
            matrix[0][j] = args.gap * j
            backtrace_matrix[0][j] = "L"
            arrow = plt.annotate(text="", xy=(j-0.2, 0), xytext=(j-0.8, 0), arrowprops=dict(arrowstyle="<-", color="r"))
            arrow_matrix[0][j]["L"] = arrow

        ax.text(j, 0, str(matrix[0][j]), va='center', ha='center')
        highlight_cell(j, 0, color="black", linewidth=1)

    for i in range(1, len(seq1)):
        if(args.algorithm == "nw"):
            matrix[i][0] = args.gap * i
            backtrace_matrix[i][0] = "U"
            arrow = plt.annotate(text="", xy=(0, i-0.15), xytext=(0, i-0.85), arrowprops=dict(arrowstyle="<-", color="r"))
            arrow_matrix[i][0]["U"] = arrow

        ax.text(0, i, str(matrix[i][0]), va='center', ha='center')
        highlight_cell(0, i, color="black", linewidth=1)

    for i in range(1, len(seq1)):
        for j in range(1, len(seq2)):
            diagonal_value = matrix[i-1][j-1] + scores[seq1[i]][seq2[j]]
            left_value = matrix[i][j-1] + args.gap
            up_value = matrix[i-1][j] + args.gap
            if(args.function == "sim"):
                matrix[i][j] = max(diagonal_value, max(left_value, up_value))
                if(args.algorithm == "sw"):
                    matrix[i][j] = max(0, matrix[i][j])
            elif(args.function == "dis"):
                matrix[i][j] = min(diagonal_value, min(left_value, up_value))

            if(args.algorithm == "nw" or (args.algorithm == "sw" and matrix[i][j] > 0)):
                if(matrix[i][j] == diagonal_value):
                    backtrace_matrix[i][j] += "D"
                    arrow = plt.annotate(text="", xy=(j-0.2, i-0.15), xytext=(j-0.8, i-0.85), arrowprops=dict(arrowstyle="<-", color="r"))
                    arrow_matrix[i][j]["D"] = arrow

                if(matrix[i][j] == left_value):
                    backtrace_matrix[i][j] += "L"
                    arrow = plt.annotate(text="", xy=(j-0.2, i), xytext=(j-0.8, i), arrowprops=dict(arrowstyle="<-", color="r"))
                    arrow_matrix[i][j]["L"] = arrow

                if(matrix[i][j] == up_value):
                    backtrace_matrix[i][j] += "U"
                    arrow = plt.annotate(text="", xy=(j, i-0.15), xytext=(j, i-0.85), arrowprops=dict(arrowstyle="<-", color="r"))
                    arrow_matrix[i][j]["U"] = arrow

                if(matrix[i][j] >= maximum_value):
                    if(matrix[i][j] > maximum_value):
                        del maximum_indices[:]

                    maximum_value = matrix[i][j]
                    maximum_indices.append([i, j])

            ax.text(j, i, str(matrix[i][j]), va='center', ha='center')
            highlight_cell(j, i, color="black", linewidth=1)

    ax.matshow(matrix, aspect='auto', cmap=ListedColormap(["w"]))
    ax.set_yticklabels([""]+list(seq1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.set_xticklabels([""]+list(seq2))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    alignments = []


    def backtrace(i, j, aligned_seq1, aligned_seq2):
        if((i == 0 and j == 0) or (args.algorithm == "sw" and matrix[i][j] == 0)):
            if(not [aligned_seq1, aligned_seq2] in alignments):
                alignments.append([aligned_seq1, aligned_seq2])

            return

        if("D" in backtrace_matrix[i][j]):
            try:
                arrow_matrix[i][j]["D"].remove()
            except ValueError as ve:
                pass

            plt.annotate(text="", xy=(j-0.2, i-0.15), xytext=(j-0.8, i-0.85), arrowprops=dict(arrowstyle="<-", color="b"))
            backtrace(i-1, j-1, [seq1[i]]+aligned_seq1, [seq2[j]]+aligned_seq2)
        if("L" in backtrace_matrix[i][j]):
            try:
                arrow_matrix[i][j]["L"].remove()
            except ValueError as ve:
                pass

            plt.annotate(text="", xy=(j-0.2, i), xytext=(j-0.8, i), arrowprops=dict(arrowstyle="<-", color="b"))
            backtrace(i, j-1, ["-"]+aligned_seq1, [seq2[j]]+aligned_seq2)
        if("U" in backtrace_matrix[i][j]):
            try:
                arrow_matrix[i][j]["U"].remove()
            except ValueError as ve:
                pass

            plt.annotate(text="", xy=(j, i-0.15), xytext=(j, i-0.85), arrowprops=dict(arrowstyle="<-", color="b"))
            backtrace(i-1, j, [seq1[i]]+aligned_seq1, ["-"]+aligned_seq2)


    if(args.algorithm == "nw"):
        highlight_cell(len(seq2)-1, len(seq1)-1, color="green", linewidth=2)
        backtrace(len(seq1)-1, len(seq2)-1, [], [])
    elif(args.algorithm == "sw" and maximum_value != 0):
        for indices in maximum_indices:
            highlight_cell(indices[1], indices[0], color="green", linewidth=2)
            backtrace(indices[0], indices[1], [], [])


    aligments_str = ""
    for alignment in alignments:
        if(not len(aligments_str)):
            aligments_str = str(alignment[0]) + "\n" + str(alignment[1])
        else:
            aligments_str += "\n\n" + str(alignment[0]) + "\n" + str(alignment[1])

    plt.gcf().text(0.93, 0.5, aligments_str, verticalalignment="center", bbox=dict(boxstyle='round', facecolor='white', alpha=0.15))
    plt.savefig("backtrace.pdf", bbox_inches="tight")
except Exception as ex:
    print(traceback.print_exc())
