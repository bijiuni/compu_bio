# Lab5 RNAfold and Riboswitches

## Overview

The aim of this lab is to run RNAfold package to obtain secondary structures of folded RNA.
The inputs are nucleic sequences and possibly constraints. All sequences and riboswitch scenarios were obatined from the paper written 
by Penchovsky and	Breaker.


## Understanding the Paper

**_Q:_** Will	the	AND-1	riboswitch	cleave	itself	when	both	of	its	OBS	are	bound?

**_A:_** Yes it will. If both OBS are bound, that means both the inputs to the AND gate are **TRUE**. The output should therefore be **TRUE**. The regions resposible of self-cleavage are aligned and bound to each other. The riboswitch will cleave itself.


**_Q:_** Will	the	OR-1	riboswitch	cleave	itself	when	neither	of	its	OBS	are	bound?

**_A:_** No it won't. If neither of the OBS are bound, that means both the inputs to the OR gate are **FALSE**. The output should therefore be **FALSE**. The regions resposible of self-cleavage are not exactly bound to each other. The riboswitch will not cleave itself.

**_Q:_** What	behavior	do	we	expect	from	the	YES-1	riboswitch?

**_A:_** If the OBS is bound, that means the input of the YES gate is **TRUE**. The output should be **TRUE** as well. The riboswitch cleaves itself.

If the OBS is not bound, that means the input of the YES gate is **FALSE**. The output should be **FALSE** as well. The riboswitch does not cleave itself.


## Grabbing the Sequences

All information obtained from the sequence is summarized into the table below:

| Riboswitch name | Start OBS-1  | End OBS-1 | Start OBS-2 | End OBS-2 | Start red1 | End red1 | Start red2 | Snd red2 |Sequence|
|------|------|------|------|------|------|------|------|------|------|
|   YES-1 | 26|47|N/A|N/A|16|21|49|54|GGGCGACCCUGAUGAGCUUGAGUUUAGCUCGUCACUGUCCAGGUUCAAUCAGGCGAAACGGUGAAAGCCGUAGGUUGCCC|
|   NOT-1  | 44|66|N/A|N/A|40|43| 74|77|GGCAGGUACAUACAGCUGAUGAGUCCCAAAUAGGACGAAACGCGACACACACCACUAAACCGUGCAGUGUUUUGCGUCCUGUAUUCCACUGC|
|   AND-1 |30|45|49|64|16|23|70|77|GGGCGACCCUGAUGAGCUUGGUUUAGUAUUUACAGCUCCAUACAUGAGGUGUUAUCCCUAUGCAAGUUCGAUCAGGCGAAACGGUGAAAGCCGUAGGUUGCCCAGAGACAAU|
|OR-1 |27|46|47|66|16|26|67|77|GGGCGACCCUGAUGAGCUUGGUUGAGUAUUUACAGCUCCAUACAUGAGGUGUUCUCCCUACGCAAGUUCGAUCAGGCGAAACGGUGAAAGCCGUAGGUUGCCC|


## Routine to fold RNA

The developed routine is shown below. We first import subprocess. When calling subprocess, we simply call the RNAfold function and use
the sequences as input in the form of ASCII characters. Terminal output and outputfile are implemented seperately.

```
import subprocess
seqs = \
""">YES-1
GGGCGACCCUGAUGAGCUUGAGUUUAGCUCGUCACUGUCCAGGUUCAAUCAGGCGAAACGGUGAAAGCCGUAGGUUGCCC
>NOT-1
GGCAGGUACAUACAGCUGAUGAGUCCCAAAUAGGACGAAACGCGACACACACCACUAAACCGUGCAGUGUUUUGCGUCCUGUAUUCCACUGC
>AND-1
GGGCGACCCUGAUGAGCUUGGUUUAGUAUUUACAGCUCCAUACAUGAGGUGUUAUCCCUAUGCAAGUUCGAUCAGGCGAAACGGUGAAAGCCGUAGGUUGCCCAGAGACAAU
>OR-1
GGGCGACCCUGAUGAGCUUGGUUGAGUAUUUACAGCUCCAUACAUGAGGUGUUCUCCCUACGCAAGUUCGAUCAGGCGAAACGGUGAAAGCCGUAGGUUGCCC
>@
"""

p = subprocess.run(['RNAfold'],
                   input=bytes(seqs,'ascii'),
                  stdout=subprocess.PIPE,
                  stderr=subprocess.PIPE)

print("*** Terminal output ***")
print(p.stderr.decode())

print("*** Output file ***")
print(p.stdout.decode())
```



## The meat of the assignment

### (1) Default parameters and the effect of temperature

With default parameters (please see ipython file for detailed codes), we obtained the following configurations.


##### YES-1

<img src=\"IMG/YES-1_ss.png\" width=\"300\">


##### NOT-1

<img src=\"IMG/NOT-1_ss.png\" width=\"300\">


##### AND-1

<img src=\"IMG/AND-1_ss.png\" width=\"300\">


##### OR-1

<img src=\"IMG/OR-1_ss.png\" width=\"300\">


**Discussion:** The obtained RNA plots are identical to those obtained by Penchovsky & Breaker for all four gates. All stem-loop structues match up so far.

By changing the temperature using code like the following, the configurations of the RNA can be changed. In order to keep the outside conditions unchanged, we used default parameters in the rest of this lab.

```
p = subprocess.run(['RNAfold', '--temp=5'],
                   input=bytes(seqs,'ascii'),
                  stdout=subprocess.PIPE,
                  stderr=subprocess.PIPE)
```


### (2) OBS binding in YES-1 and NO-1

To automate the constraint generating process, and since the OBS positions have already been identified in the beginning of this lab, we use the following code as a template to simulate RNA folding in all binding conditions.

```
seq = "GGGCGACCCUGAUGAGCUUGAGUUUAGCUCGUCACUGUCCAGGUUCAAUCAGGCGAAACGGUGAAAGCCGUAGGUUGCCC"
bstart = 26
bend = 47
end = len(seq)

string = ">YES-1_T\n" +seq+"\n"

for i in range(bstart):
    string += "."    
for i in range(bstart, bend+1):
    string += "x"    
for i in range(bend+1, end):
    string += "." 

string += "\n>@"
    
p = subprocess.run(['RNAfold', '-C'],
                       input=bytes(string, 'ascii'),
                      stdout=subprocess.PIPE,
                      stderr=subprocess.PIPE) 
```

The codes are self-explanatory. The input is the sequence of the RNA followed by the '.x' constraint. With the -C parameters added, the RNAfold function will generate the RNA folding configurations in the constraint conditions.



Here are the results:


#### YES-1 constrained folding

>YES-1_T
GGGCGACCCUGAUGAGCUUGAGUUUAGCUCGUCACUGUCCAGGUUCAAUCAGGCGAAACGGUGAAAGCCGUAGGUUGCCC
((((((((.......((((((...........................))))))...(((((....))))).)))))))) (-24.50)

<img src=\"IMG/YES-1_T_ss.png\" width=\"300\">


#### NOT-1 constrained folding

>NOT-1_T
GGCAGGUACAUACAGCUGAUGAGUCCCAAAUAGGACGAAACGCGACACACACCACUAAACCGUGCAGUGUUUUGCGUCCUGUAUUCCACUGC
.((((....((((((..((((.((((......))))(((((...........................))))).))))))))))....)))) (-16.00)

<img src=\"IMG/NOT-1_T_ss.png\" width=\"300\">


**Discussion:** For yes-1 riboswitch, the RNA configurations are exactly the same with and without constraint. If the OBS is bound, that means the input of the YES gate is **TRUE**. The output is **TRUE** as well. The riboswitch cleaves itself. If the OBS is not bound, that means the input of the YES gate is **FALSE**. The output is **FALSE** as well. The riboswitch does not cleave itself. We therefore can generate the truth table as below.


#### YES-1 truth table


|Input|Output|
|------|------|
|0|0|
|1|1|


for not-1 riboswitch, the RNA folding result without constraint is the same with the paper as previously showed. However, when it comes to constraint case, a pair of GU and a pair of CG are missing at position 42-43, 67-68. Despite this, the two red regions are still not aligned, and the output is still **FALSE**. We therefore can generate the truth table as below.


#### NOT-1 truth table


|Input|Output|
|------|------|
|0|1|
|1|0|


### (3) OBS binding in AND-1 and OR-1

Using the same template and parameters (detailed codes please see the ipython file), we generated the RNA folding of AND-1 and OR-1 in T/F, F/T, and T/T situations respectively. Here are the results:

#### AND-1 OBS-1 bound alone

>AND-1_TF
GGGCGACCCUGAUGAGCUUGGUUUAGUAUUUACAGCUCCAUACAUGAGGUGUUAUCCCUAUGCAAGUUCGAUCAGGCGAAACGGUGAAAGCCGUAGGUUGCCCAGAGACAAU
((((((((((((((((((((..........................(((.(....))))...))))))).))))).....(((((....))))).))))))))......... (-33.90)

<img src=\"IMG/AND-1_TF_ss.png\" width=\"300\">


##### AND-1 OBS-2 bound alone

>AND-1_FT
GGGCGACCCUGAUGAGCUUGGUUUAGUAUUUACAGCUCCAUACAUGAGGUGUUAUCCCUAUGCAAGUUCGAUCAGGCGAAACGGUGAAAGCCGUAGGUUGCCCAGAGACAAU
((((((((((.(((....(((...(((.......))))))..))).))..................((((......))))(((((....))))).))))))))......... (-28.30)

<img src=\"IMG/AND-1_FT_ss.png\" width=\"300\">


##### AND-1 OBS-1 and OBS-2 both bound

>AND-1_TT
GGGCGACCCUGAUGAGCUUGGUUUAGUAUUUACAGCUCCAUACAUGAGGUGUUAUCCCUAUGCAAGUUCGAUCAGGCGAAACGGUGAAAGCCGUAGGUUGCCCAGAGACAAU
(((((((((((((((((................................................)))).))))).....(((((....))))).))))))))......... (-26.30)

<img src=\"IMG/AND-1_TT_ss.png\" width=\"300\">

**Discussion:** For and-1 riboswitch, four situations are discussed:
FF: Exactly the same with the paper. Both blue areas are not bound, the red areas are not aligned, the riboswitch does not undergo cleavage.

TF: The configurations are basically the same with the paper, except that in our case position 28-56 are a few stem loops instead of a large stem loop in the paper. The two red regions are not aligned. The ouput is false. The riboswitch does not undergo cleavage.

FT: This case is a bit messed up. A few stem loops are not the same with the paper. For example, at position 10-40 and 65-80. However, the two red regions are still not aligned. The result is still false.

TT: Exactly the same with the paper. Both blue areas are bound, the red areas are aligned and bound, the riboswitch will undergo cleavage. The ouput is true in this case.

We therefore can generate the truth table as below.


#### AND-1 truth table

|Input 1|Input 2|Output|
|------|------|------|
|0|0|0|
|1|0|0|
|0|1|0|
|1|1|1|




#### OR-1 OBS-1 bound alone

>OR-1_TF
GGGCGACCCUGAUGAGCUUGGUUGAGUAUUUACAGCUCCAUACAUGAGGUGUUCUCCCUACGCAAGUUCGAUCAGGCGAAACGGUGAAAGCCGUAGGUUGCCC
((((((((((((((((((((((.(((...........................)))...)).))))))).))))).....(((((....))))).)))))))) (-34.20)

<img src=\"IMG/OR-1_TF_ss.png\" width=\"300\">


##### OR-1 OBS-2 bound alone

>OR-1_FT
GGGCGACCCUGAUGAGCUUGGUUGAGUAUUUACAGCUCCAUACAUGAGGUGUUCUCCCUACGCAAGUUCGAUCAGGCGAAACGGUGAAAGCCGUAGGUUGCCC
((((((((.......((((((((((..........................................))))))))))...(((((....))))).)))))))) (-28.84)

<img src=\"IMG/OR-1_FT_ss.png\" width=\"300\">


##### OR-1 OBS-1 and OBS-2 both bound

>OR-1_TT
GGGCGACCCUGAUGAGCUUGGUUGAGUAUUUACAGCUCCAUACAUGAGGUGUUCUCCCUACGCAAGUUCGAUCAGGCGAAACGGUGAAAGCCGUAGGUUGCCC
((((((((.......((((((((((..........................................))))))))))...(((((....))))).)))))))) (-28.84)

<img src=\"IMG/OR-1_TT_ss.png\" width=\"300\">

**Discussion:** For OR-1 riboswitch, four situations are discussed:
FF: Exactly the same with the paper. Both blue areas are not bound, the red areas are not aligned, the riboswitch does not undergo cleavage.

TF: This case is a bit messed up. A few stem loops are not the same with the paper. In addition, the two red regions are not aligned at all. The result is still false.

FT: In terms of the red region, most of the base pairs are aligned and bound expect the GU pair at position 67 and 26. Still we count this output as false.

TT: Not exactly the same with the paper. A few stem loops are messed up. However, when both blue areas are bound, the red areas are aligned and bound, the riboswitch will undergo cleavage. The ouput is true in this case.

We therefore can generate the truth table as below.

#### OR-1 truth table

|Input 1|Input 2|Output|
|------|------|------|
|0|0|0|
|1|0|0|
|0|1|0|
|1|1|1|




## Final Question: According to your results, do the riboswitches work as the paper claims?

If we neglect the different configurations and only focus on the red regions' alignment, three riboswitches work as the paper claims: YES-1, NOT-1, and AND-1. The TF and FT situations in OR-1 gate does not work as claimed.

Both the change in configurations and the different output in OR-1 riboswitch can probabaly be explained by the different parameters in the RNAfold function. As can be seen in the document, a number of parameters including temperature can be adjusted. They can all lead to different RNA folding configurations. One thing can be learned here is that specifications of technical details like parameters are important if one wants his or her research to be repeatable.
