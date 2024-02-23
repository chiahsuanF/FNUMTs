#!/usr/bin/env python3

import pandas as pd

def merge_intervals(bed_list, max_distance):
    if not bed_list:
        return []

    bed_list = sorted(bed_list, key=lambda x: x[1])
    merged = []
    current_interval = bed_list[0]
    count = 1  

    for i in range(1, len(bed_list)):
        if bed_list[i][0] == current_interval[0] and bed_list[i][1] <= current_interval[2] + max_distance:
            current_interval = (current_interval[0], current_interval[1], max(current_interval[2], bed_list[i][2]))
            count += 1
        else:
            merged.append((current_interval[0], current_interval[1], current_interval[2], count))
            current_interval = bed_list[i]
            count = 1
    merged.append((current_interval[0], current_interval[1], current_interval[2], count))
    return merged

def read_intervals(filename):
    intervals = pd.read_csv(filename, sep='\t', header=None).values
    return sorted(intervals, key=lambda x: (x[0], x[1], x[2]))

def bed_intersect(a_intervals, b_intervals):
    intersecting_intervals = []
    i, j = 0, 0

    while i < len(a_intervals) and j < len(b_intervals):
        a_interval = a_intervals[i]
        b_interval = b_intervals[j]
        
        if a_interval[0] < b_interval[0] :
            i += 1
        elif a_interval[0] > b_interval[0] :
            j += 1       
        elif (a_interval[0] == b_interval[0])  and  (a_interval[1] < b_interval[2]) and ((a_interval[2]) >= (b_interval[1])):
            intersecting_intervals.append(list(a_interval) + list(b_interval))
            j += 1             
        elif (a_interval[0] == b_interval[0]) and ((a_interval[2] ) < (b_interval[1])):
            i += 1
        else:
            j += 1
            
    return intersecting_intervals
