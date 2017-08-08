# LeeStrazicichUnitRoot
R code to perform the Lee Strazicich Unit Root test by Lee, Strazicich (2003, 2004), which allows
the endogenous determination of one or two structural breaks.
The code is based on the original code by Junsoo Lee and Mark C. Strazicich additionally I took the available RATS code [1] by Tom Doan and replicated my results with this procedure.
I also now implemented an updated version of my code, which allows to run the most computationally demanding functions in parallel. This allows a really fast calculation of the results, compared with the first implementation in R. Additionally I published the procedures to calculate bootstrapped results for the test statistic and the break dates. This is mostly based on the paper by Chou (2007).
For an example of the usage and the presentation of this test, based on the code, please be referred to the published paper or the previous conference paper version. [2,3]

Usage:

1) The R files LeeStrazizichUnitRootTest.R and LeeStrazizichUnitRootTestParallelization.R contain the functions, which need to be sourced in the R environment.

2) The R file LeeStrazizichApplication.R then just gives two quick examples on how to use the two functions and get the results.


I would like to thank the authors, who made the code available to me.

[1] https://estima.com/forum/viewtopic.php?f=7&t=126

Lee, Junsoo and Mark C. Strazicich (2003). "Minimum Lagrange Multiplier Unit
Root Test with Two Structural Breaks". In: The Review of Economics and Statistics 85.4, pp. 1082-1089.

Lee, Junsoo and Mark C. Strazicich (2004). "Minimum LM Unit Root Test with One Structural Break". In: 04-17. url: https://ideas.repec.org/p/apl/wpaper/04-17.html (visited on 02/04/2015).

Chou, Win Lin (2007). "Performance of LM-type unit root tests with trend break: A bootstrap approach", Economics Letters, Volume 94, Issue 1, January 2007, Pages 76-82, ISSN 0165-1765, http://dx.doi.org/10.1016/j.econlet.2006.08.004.
(http://www.sciencedirect.com/science/article/pii/S0165176506002722)

[2] Lips, J. (2017). Do They Still Matter? – Impact of Fossil Fuels on Electricity Prices in the Light of Increased Renewable Generation. Journal of Time Series Econometrics, 9(2), pp. -. Retrieved 8 Aug. 2017 from https://doi.org/10.1515/jtse-2016-0018

[3] Lips, Johannes (2016). "Do They Still Matter? - Impact of Fossil Fuels on Electricity Prices in the Light of Increased Renewable Generation." In . Beiträge Zur Jahrestagung Des Vereins Für Socialpolitik 2016: Demographischer Wandel - Session: Auctions and Prices. Kiel und Hamburg: ZBW - Deutsche Zentralbibliothek für Wirtschaftswissenschaften, Leibniz-Informationszentrum Wirtschaft. http://hdl.handle.net/10419/145601.
