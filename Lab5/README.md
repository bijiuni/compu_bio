# Lab5 RNAfold and Riboswitches

## Overview

The aim of this lab is to run RNAfold package to obtain secondary structures of folded RNA.
The inputs are nucleic sequences and possibly constraints. All sequences and riboswitch scenarios were obatined from the paper written 
by Penchovsky and	Breaker.

## Understanding the Paper

*Q* Will	the	AND-1	riboswitch	cleave	itself	when	both	of	its	OBS	are	bound?
Yes it will.


Will	the	OR-1	riboswitch	cleave	itself	when	neither	of	its	OBS	are	bound?


What	behavior	do	we	expect	from	the	YES-1	riboswitch?


## Grabbing the Sequences

All information is summarized into the table below:



YES-1: (Figure 2a)
```
seq:GGGCGACCCUGAUGAGCUUGAGUUUAGCUCGUCACUGUCCAGGUUCAAUCAGGCGAAACGGUGAAAGCCGUAGGUUGCCC
OBS: 26-47
Red regions: 16-21, 49-54
```

NOT-1: (Figure 4a)
```
seq:GGCAGGUACAUACAGCUGAUGAGUCCCAAAUAGGACGAAACGCGACACACACCACUAAACCGUGCAGUGUUUUGCGUCCUGUAUUCCACUGC
OBS: 44-66
Red regions: 40-43, 74-77
```

AND-1: (Figure 5a)
```
seq:GGGCGACCCUGAUGAGCUUGGUUUAGUAUUUACAGCUCCAUACAUGAGGUGUUAUCCCUAUGCAAGUUCGAUCAGGCGAAACGGUGAAAGCCGUAGGUUGCCCAGAGACAAU
OBS1: 30-45
OBS2: 49-64
Red regions: 16-23, 70-77
```

OR-1: (Figure 6a)
```
seq:GGGCGACCCCUGAUGGCUUGGUUGAGUAUUUACAGCUCCAUAUACAUGAGGUGUUCUCCCUACGCAAGUUCGAUCAGGCGAAACGGUGAAAGCCGUAGGUUGCCC
OBS1: 26-4 
OBS2: 47-66
Red regions: 16-26, 67-77
```
|Riboswitch|Sequence|Start	and	end	coordinates	of	OBS-1 (blue region)|Start and end coordinates of OBS-2 (blue	region, only applicable	to AND-1 and OR-1)|	Start and end coordinates of the two red regions|
|:----------|:--------|:----------|:--------|:--------|
|YES-1|GGGCGACCCUGAUGAGCUUGAGUUUAGCUCGUCACUGUCCAGGUUCAAUCAGGCGAAACGGUGAAAGCCGUAGGUUGCCC|26-47| None |16-21, 49-54|
|NOT-1|GGCAGGUACAUACAGCUGAUGAGUCCCAAAUAGGACGAAACGCGACACACACCACUAAACCGUGCAGUGUUUUGCGUCCUGUAUUCCACUGC|44-66| None | 40-43, 74-77|
|AND-1|GGGCGACCCUGAUGAGCUUGGUUUAGUAUUUACAGCUCCAUACAUGAGGUGUUAUCCCUAUGCAAGUUCGAUCAGGCGAAACGGUGAAAGCCGUAGGUUGCCCAGAGACAAU|30-45|49-64|16-23, 70-77|
|OR-1|GGGCGACCCCUGAUGGCUUGGUUGAGUAUUUACAGCUCCAUAUACAUGAGGUGUUCUCCCUACGCAAGUUCGAUCAGGCGAAACGGUGAAAGCCGUAGGUUGCCC|27-46|47-66|16-26, 67-77|
