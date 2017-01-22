import csv
import numpy as np


def ReadColumn(fileName):
    points = []
    with open(fileName, 'rb') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',')
       	for row in csvreader:
       		points.append(map(int, row))
    return np.array(points)
