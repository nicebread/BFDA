\name{BFDA}
\title{News for package 'BFDA'}
\encoding{UTF-8}

\section{Version 0.5}{

Thanks to Tom Heyman for many helpful suggestions!

\subsection{New features}{
	\itemize{
		\item Added BFDA for Bayesian A/B tests (see argument type = "abtest" in the BFDA.sim() function)
		\item Plot for fixed-N design can now be drawn separately for H0 and H1 (see evDens() function: only one out of the arguments BFDA.H1 and BFDA.H0 has to be defined)
	}
}

\subsection{Bug fixes/Glitches}{
	\itemize{
		\item Fixed bug in summary print of an BFDA.sim object: Mean and quantiles of simulated (random) true effect sizes now are correctly displayed.
	}
}

}


