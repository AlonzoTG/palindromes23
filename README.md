# palindromes23
Some useless number theory code to find the elusive dual palindromes in Binary and Ternary. 

Most of this code was done in late 2013 and early 2014.

This code finds numbers that are palindromes in BOTH binary and ternary. This set would be entirely ho-hum if it weren't so damn sparse! The distribution is fairly random, there is some controvercy as to whether there is any pattern in the distribution. However, there are roughly 4 orders of magnitude difference between members of this set (the distribution is exponential!), only 17 examples have been found, half of which were fonud on my previous computer. 

Some other ppl came around later and claimed some of my numbers, BLEH!!! I found them first and I have proof! 

REQUIREMENTS IN CURRENT CONFIGURATION. 

I am running this code right now, looking at new search space. Therefore I initiated the main variable "i" in main2.c with the edge of the search space. I expect it to take my Ryzen R7 1800x about six months to complete this block. If you want to study how the program works or test it, set it to 1 and see the list from the beginning. 

This code requires an AUTHENTICAMD processor. I use a few instructions that are not implemented on intel chips, last I checked. 
It is known to run on linux.
It is currently configured to run on 32GB of ram. (uses about 25gb) 
The GMP library is required, I did customize one of the headers (provided) but it still seems to work fine. 
It is currently configured to use 15 threads of execution, I suggest configuring it for one less than the number of hardware threads. 


Tuning.

In order to tune this code you need to understand how it works. 

Lets say we have a ternary palindrome of the form: 

ABCDEXXXXXXXXYYYYYYY1YYYYYYYXXXXXXXEDCBA

We only spend time thinking about the top half because the bottom half is always the same.... 

So our counter reads: 

ABCDEXXXXXXX  

We know the end we want will be EDCBA. 

We also know that with XXXXXXX, our end will be ZZZZZ. 

So we subtract ZZZZZ from EDCBA and look that up in a great big 24gb table and get a list of not more than 11 entries we have to look at, so we go from 500,000,000 numbers to a hotlist of just 11 (!!!) "Eviltable.txt" is a sample of this data...

Now 999,999,999,999 times out of a trillion, the XXXXXXX will turn out to be wrong, try again, but still the speedup is real! 

If we could go to, say, 100gb of ram, we could quadruple the speedup! 

The search through each of the blocks is embarasingly paralell, so go ahead and unleash a threadripper or an EPYC on this code... The dispatch code, in the critical section, is a bit dodgy though and needs to be verified a few more times. =\ 

So you need to configure the table size, I'm almost maxed out before I have to rewrite stuff (!!), and then ajust the number of ternary digits you are masking, (I forget exactly how it works, not sure the existing code is correct). When you run it, it will report the maximum entry size in the table, adjust the width of the table to that size for the best ram usage, currently 11. 

My favorite make targets: 

make debug  -> single threaded debug version.
make stable -> let'er rip! 
make clean  -> delete all detritus. 


TODO: 

-> Make the "squelch" variable adaptive to the value of i
-> verify the code doesn't have any blind spots. 
-> any further optimizations. 
