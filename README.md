bwa-spark-0.2.0
===============
Development NOTE:
(1) The order after sorting INFLUENCES the results. The result will be slightly different from the original C version.
    This occurs in MemChainToAlign(), MemSortAndDedup() and MemMarkPrimarySe().
(2) Tree ordering influences results... (B-tree: original C implementation  V.S. RB-tree: our implementation)
