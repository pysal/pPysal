"""
merge sort parallelization
"""

import math


def mergeSort(a):
    len_a = len(a)
    if len_a <= 1:
        return a
    m = int(math.floor(len_a) / 2)
    left = a[0:m]
    right = a[m:]
    left = mergeSort(left)
    right = mergeSort(right)
    return merge(left, right)

def merge(left, right):
    a = []
    while len(left) > 0 or len(right) > 0:
        if len(left) > 0 and len(right) > 0:
            if left[0] <= right[0]:
                a.append(left.pop(0))
            else:
                a.append(right.pop(0))
        elif len(left) > 0:
            a.append(left.pop(0))
        elif len(right)> 0:
            a.append(right.pop(0))
    return a


if __name__ == '__main__':

    l = [ 1, 2, 7]
    r = [ 3, 5, 20]

    s = l[:]
    s.extend(r[:])

    print merge(l, r)
    print s
    print mergeSort(s)
