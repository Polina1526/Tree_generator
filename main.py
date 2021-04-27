import argparse
import functions as fun
from tree import Tree


parser = argparse.ArgumentParser()
parser.add_argument("-ac", "--average_calculating", type=int, default=1,
                    help="generate '-ac' different trees and calculate the average TMRCA")
parser.add_argument("-txt", "--txt_file", type=bool, default=False,
                    help="whether it is necessary to output history to a txt file"
                         "with name 'tree_history' (default=False)")
parser.add_argument("-d", "--draw_tree", type=bool, default=False,
                    help="will the tree be drawn (default=False),"
                         "when drawing a tree, it is better to use T = -1, since otherwise the result is undefined")
args = parser.parse_args()


tree = Tree("coefficients", "migration_rates")


tmrca = []
for i in range(args.average_calculating):
    tmrca.append(fun.generate_tree(tree, args.txt_file))
    if args.draw_tree:
        fun.draw_tree(tree)

if args.average_calculating > 1:
    tmrca_average = sum(tmrca) / args.average_calculating
    print(tmrca_average)
