# LeeStrazicichUnitRoot
R code to perform the Lee Strazicich Unit Root test by Lee, Strazicich (2003, 2004), which allows
the endogenous determination of one or two structural breaks.
The code is based on the original code by Junsoo Lee and Mark C. Strazicich additionally I took the available RATS code [1] by Tom Doan and replicated my results with this procedure.
I also now implemented an updated version of my code, which allows to run the most computationally demanding functions in parallel. This allows a really fast calculation of the results, compared with the first implementation in R. Additionally I am happy to provide code to produce bootstrapped confidence intervals for the test statistic and the break dates. The code is not in publishable form and therefore not yet uploaded to github. This is based on the paper by Chou (2007).

Usage:

1) The R files LeeStrazizichUnitRootTest.R and LeeStrazizichUnitRootTestParallelization.R contain the functions, which need to be sourced in the R environment.

2) The R file LeeStrazizichApplication.R then just gives two quick examples on how to use the two functions and get the results.


I would like to thank the authors, who made the code available to me.

[1] https://estima.com/forum/viewtopic.php?f=7&t=126

Lee, Junsoo and Mark C. Strazicich (2003). "Minimum Lagrange Multiplier Unit
Root Test with Two Structural Breaks". In: The Review of Economics and Statistics 85.4, pp. 1082-1089.

Lee, Junsoo and Mark C. Strazicich (2004). "Minimum LM Unit Root Test with One Structural Break". In: 04-17. url: https://ideas.repec.org/p/apl/wpaper/04-17.html (visited on 02/04/2015).

Win Lin Chou (2007). "Performance of LM-type unit root tests with trend break: A bootstrap approach", Economics Letters, Volume 94, Issue 1, January 2007, Pages 76-82, ISSN 0165-1765, http://dx.doi.org/10.1016/j.econlet.2006.08.004.
(http://www.sciencedirect.com/science/article/pii/S0165176506002722)
