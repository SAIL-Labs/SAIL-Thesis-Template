%---------------------------------------------------------------
% PhDthesis version 2 style file by Jakob Suckale, 2007
% based on CUEDthesis version 1 by Harish Bhanderi 2002
%---------------------------------------------------------------


%-------------------------- identification ---------------------
%\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{Latex/Classes/PhDthesisPSnPDF}[2007/09/06 v2 PhD thesis class]

%:-------------------------- report or book -----------------------

%If you want to use a Report style document then uncomment the following 3 lines and comment the below 8 book style lines:
%\DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}
%\ProcessOptions\relax
%%\LoadClass[a4paper]{report}
%\ifx\pdfoutput\undefined
%   \LoadClass[dvips, a4paper]{report}
%\else
%   \LoadClass[pdftex, a4paper]{report}
%\fi

%%Bibliography
%%uncomment next line to change bibliography name to references for Report document class
%\renewcommand{\refname}{References}

%If you want to use a Book style document then uncomment the following 3 lines and comment the above 8 report style lines:
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{book}}
\ProcessOptions\relax

\LoadClass[pdftex, a4paper]{book}


%:-------------------------- packages for fancy things -----------------------
\newif\ifthesisisdraft
\usepackage{amssymb}
\usepackage{amsmath}
%\usepackage{wrapfig} % to include figure with text wrapping around it
\usepackage[margin=10pt,font=small,labelfont=bf]{caption} % for improved layout of figure captions with extra margin, smaller font than text
\usepackage{fancyhdr} % for better header layout
\usepackage{eucal}
\usepackage[english]{babel}
\usepackage[usenames, dvipsnames]{color}
\usepackage[perpage]{footmisc}
%\usepackage[square, sort, numbers, authoryear]{natbib}
\usepackage[autostyle]{csquotes}
\usepackage{ifthen}
\usepackage{multicol
} % for pages with multiple text columns, e.g. References
\setlength{\columnsep}{20pt} % space between columns; default 10pt quite narrow
%\setlength{\multicolsep}{6.0pt plus 2.0pt minus 1.5pt}

\usepackage[nottoc]{tocbibind} % correct page numbers for bib in TOC, nottoc suppresses an entry for TOC itself
%\usepackage{nextpage}
\usepackage{paralist}
\usepackage{afterpage}
\usepackage{xspace}
\usepackage{booktabs}
\usepackage{gensymb}
\usepackage{subcaption}

%:-------------------------- PDF/PS setup -----------------------

%Assuming always PDF output, no longer need ifpdf package
%\usepackage{ifpdf} 

%if you use an older version, uncomment these lines:
%\newif \ifpdf
%    \ifx \pdfoutput \undefined
        % for running latex
%        \pdffalse
%    \else
        % for running pdflatex
%        \pdfoutput = 1    % positive value for a PDF output
                          % otherwise a DVI
%        \pdftrue
%\fi

% PDF output configuration (assuming always PDF)
%-->
%--> Google.com search "hyperref options"
%--> 
%--> http://www.ai.mit.edu/lab/sysadmin/latex/documentation/latex/hyperref/manual.pdf
%--> http://www.chemie.unibas.ch/~vogtp/LaTeX2PDFLaTeX.pdf 
%--> http://www.uni-giessen.de/partosch/eurotex99/ oberdiek/print/sli4a4col.pdf
%--> http://me.in-berlin.de/~miwie/tex-refs/html/latex-packages.html
%-->
% \PassOptionsToPackage{hyphens}{url}
\usepackage[ pdftex, plainpages = false, pdfpagelabels, 
             pdfpagelayout = OneColumn, % display single page, advancing flips the page
             bookmarks,
             bookmarksopen = true,
             bookmarksnumbered = true,
             breaklinks = true,
             linktocpage,
             pagebackref=false,
             colorlinks = true,  % was false % change for printing. #todo
             linkcolor = blue,
             urlcolor  = blue,
             citecolor = red,
             anchorcolor = green,
             hyperindex = true,
             hyperfigures,
             backref=false,
             ]{hyperref} 

\ifthesisisdraft
    \usepackage[pdftex,draft]{graphicx}
\else
    \usepackage[pdftex]{graphicx}
\fi
\DeclareGraphicsExtensions{.png, .jpg, .jpeg, .pdf, .gif} %GIF doesn't work??
\pdfcompresslevel=9
\graphicspath{{0_frontmatter/figures/PNG/}{0_frontmatter/figures/PDF/}{0_frontmatter/figures/}}

%: --------------------------------------------------------------
%:               Glossary/Abbrev./Symbols 
%: --------------------------------------------------------------
\usepackage[toc,nonumberlist]{glossaries}
\usepackage{glossary-mcols}


%: --------------------------------------------------------------
%:                     biblatex
%: --------------------------------------------------------------
% todo biblatex entry
\usepackage[backref=true,
           hyperref=true,
            backend=biber,
              style=authoryear-comp,
            sorting=nyt,
         uniquename=false,
                url=false,
        maxbibnames=99,
       maxcitenames=2,
         uniquelist=false,
         giveninits=false,%refsegment=chapter
         ]{biblatex}
            
\defbibheading{subbibliography}{\section*{Chapter \ref{refsegment:\therefsection\therefsegment}}}

\renewbibmacro{in:}{%
    \ifentrytype{article}{}{\printtext{\bibstring{in}\intitlepunct}}}
    

%\makeatletter
\newcommand{\tempmaxup}[1]{\def\blx@maxcitenames{\blx@maxbibnames}#1}
%\makeatother

\DeclareCiteCommand{\longfullcite}[\tempmaxup]
  {\usebibmacro{prenote}}
  {\usedriver
     {\DeclareNameAlias{sortname}{default}}
     {\thefield{entrytype}}}
  {\multicitedelim}
  {\usebibmacro{postnote}}

\DeclareBibliographyCategory{fullcited}
\newcommand{\fullciteandexclude}[1]{\longfullcite{#1}\addtocategory{fullcited}{#1}}


% possessive citation - simple fallback for modern biblatex compatibility
% Define simple fallback commands that work with modern biblatex
\newrobustcmd*{\posscite}[1]{\textcite{#1}'s}
\newrobustcmd*{\Posscite}[1]{\Textcite{#1}'s}
\newrobustcmd*{\posscites}[1]{\textcites{#1}'s}


%: --------------------------------------------------------------
%:                    other defs
%: --------------------------------------------------------------

% Stolen from cwebmac.tex
% \datethis to get the date and time on the output.
%
% modified 3/27/92 (sanjeev@cs.psu.edu)
%

\def\today{\ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space\number\day, \number\year}

\newcount\twodigits
\newcount\hrs

\def\hours{\twodigits=\time%
  \ifnum\twodigits>719\gdef\ampm{pm}%
  \else\gdef\ampm{am}\fi%
  \divide\twodigits by 60%
  \hrs=\twodigits%
  \ifnum\twodigits=0\twodigits=12\fi%
  \ifnum\twodigits>12\advance\twodigits-12\fi\advance\twodigits0%
  \ifnum\twodigits>9\printtwodigits:\else\printonedigit:\fi%
  \twodigits=\hrs%
  \multiply\twodigits by-60%
  \advance\twodigits by\time%
  \printtwodigits\space\ampm}

\def\gobbleone1{}
\def\gobbleten10{}

\def\printonedigit{%
  \advance\twodigits100
  \expandafter\gobbleten\number\twodigits
  \advance\twodigits-100 }

\def\printtwodigits{%
  \advance\twodigits100
  \expandafter\gobbleone\number\twodigits
  \advance\twodigits-100 }

\def\datethis{\medskip{\sc\today\ at \hours}\bigskip}

%:-------------------------- page layout -----------------------

%A4 settings
% A4 settings for PDF
\pdfpageheight=297mm
\pdfpagewidth=210mm

\setlength{\hoffset}{0.00cm}
\setlength{\voffset}{0.00cm}

%: Uncomment this section for two-sided printing
% ------------------------------
%\setlength{\oddsidemargin}{1.5cm}
%\setlength{\evensidemargin}{0cm}
%\setlength{\topmargin}{1mm}
%\setlength{\headheight}{1.36cm}
%\setlength{\headsep}{1.00cm}
%\setlength{\textheight}{20.84cm}
%\setlength{\textwidth}{14.5cm}
%\setlength{\marginparsep}{1mm}
%\setlength{\marginparwidth}{3cm}
%\setlength{\footskip}{2.36cm}


%: Uncomment this section for one-sided printing
% taken from the original file, but with the first two lanes modified
% ------------------------------
%\setlength{\evensidemargin}{1.9cm} % was 1.96cm in original
%\setlength{\oddsidemargin}{-0.001cm} % was -0.54cm in original file
%\setlength{\topmargin}{1mm}
%\setlength{\headheight}{1.36cm}
%\setlength{\headsep}{1.00cm}
%\setlength{\textheight}{20.84cm}
%\setlength{\textwidth}{14.5cm}
%\setlength{\marginparsep}{1mm}
%\setlength{\marginparwidth}{3cm}
%\setlength{\footskip}{2.36cm}


%% for digital
\setlength{\evensidemargin}{1.9cm} % was 1.96cm in original
\setlength{\oddsidemargin}{0.71cm} % was -0.54cm in original file
\setlength{\topmargin}{1mm}
\setlength{\headheight}{1.36cm}
\setlength{\headsep}{1.00cm}
\setlength{\textheight}{20.84cm}
\setlength{\textwidth}{14.5cm}
\setlength{\marginparsep}{1mm}
\setlength{\marginparwidth}{3cm}



%\RequirePackage[textwidth=14.5cm, 
%			   textheight=20.84cm,
%			   left=3.25cm,
%			   a4paper]{geometry}



%: section below defines fancy page layout options
% ------------------------------
\pagestyle{fancy}
\renewcommand{\chaptermark}[1]{\markboth{\MakeUppercase{\thechapter. #1 }}{}}
\renewcommand{\sectionmark}[1]{\markright{\thesection\ #1}}
\fancyhf{}
\fancyhead[RO]{\bfseries\rightmark}
\fancyhead[LE]{\bfseries\leftmark}
\fancyfoot[C]{\thepage}
\renewcommand{\headrulewidth}{0.5pt}
\renewcommand{\footrulewidth}{0pt}
\addtolength{\headheight}{0.5pt}
\fancypagestyle{plain}{
  \fancyhead{}
  \renewcommand{\headrulewidth}{0pt}
}



%:-------------------------- title page layout -----------------------

% starts roman page numbering until chapter 1
% important to avoid two pages numbered 1 and 2 which may cause bad links
% bug: cover i + back side ii and then numbering restarts with i; should be iii
\renewcommand{\thepage}{\roman{page}}

\newcommand{\submittedtext}{{A thesis submitted for the degree of}}

% DECLARATIONS
% These macros are used to declare arguments needed for the
% construction of the title page and other preamble.

% The year and term the degree will be officially conferred
\def\degreedate#1{\gdef\@degreedate{#1}}
\def\degreeyear#1{\gdef\@degreeyear{#1}}

% The full (unabbreviated) name of the degree
\def\degree#1{\gdef\@degree{#1}}
% The name of your college or department(eg. Trinity, Pembroke, Maths, Physics)
\def\collegeordept#1{\gdef\@collegeordept{#1}}
% The name of your University
\def\university#1{\gdef\@university{#1}}
% Defining the crest
\def\crest#1{\gdef\@crest{#1}}

% Accessor commands with error checking
\newcommand{\thetitle}{\@title}
\newcommand{\theauthor}{\@author}
\newcommand{\thedate}{\@date}

% These macros define an environment for front matter that is always 
% single column even in a double-column document.

\newenvironment{alwayssingle}{%
       \@restonecolfalse\if@twocolumn\@restonecoltrue\onecolumn
       \else\newpage\fi}
       {\if@restonecol\twocolumn\else\newpage\fi}

%%define title page layout
%\renewcommand{\maketitle}{%
%\begin{alwayssingle}
%    \renewcommand{\footnotesize}{\small}
%    \renewcommand{\footnoterule}{\relax}
%    \thispagestyle{empty}
%%  \null\vfill
%  \begin{center}
%    
%%{\large  {{\@crest} \par}}
%%{\Huge \vspace*{25mm} {\bfseries {\@title}} \par \vspace*{25mm}}
%%{{\Large \@author} \par}
%
%{\large  {{\@crest} \par}}
%{\Huge \vspace*{25mm} {\bfseries \textcolor{midgrey}{\@title}} \par \vspace*{25mm}}
%{{\Large \@author} \par}
%
%%	\vspace*{1ex}
%%	{{\@cityofbirth} \par}
%{\large 
%	\vspace*{1ex}
%    {{\@collegeordept} \par}
%	\vspace*{1ex}
%    {{\@university} \par}
%	\vspace*{25mm}
%    {{\submittedtext} \par}
%	\vspace*{1ex}
%    {\it {\@degree} \par}
%	\vspace*{2ex}
%    {\@degreedate}
%}%end large
%  \end{center}
%  \null\vfill
%\end{alwayssingle}}

% page number for cover back side should have page number blanked


\renewcommand{\maketitle}{
	\thispagestyle{empty}
	\vspace*{\fill}
	\vspace{75pt}
	\begin{center}
	\Huge \textcolor{midgrey}{\@title} \normalsize \\
	\vspace{25pt}
	by \\
	\vspace{15pt}
	\@author\\
	\vspace{50pt}
	{\large  {{\@crest} \par}}
	\vspace{50pt}
	\textsc{a thesis submitted in fulfilment\\
	of the requirements for the degree of \\
	\vspace{8pt}
	\@degree \\
	\@collegeordept \\
	\@university \\
	\vspace{10pt}
	\@degreedate}
	\end{center}
	\vspace*{\fill}
}



%:-------------------------- front matter layout -----------------------

% DEDICATION
%
% The dedication environment makes sure the dedication gets its
% own page and is set out in verse format.

\newenvironment{dedication}
{\begin{alwayssingle}
  \pagestyle{empty}
  \begin{center}
  \vspace*{1.5cm}
  {\LARGE }
  \end{center}
  \vspace{0.5cm}
  \begin{quote} \begin{center}}
{\end{center} \end{quote} \end{alwayssingle}}


% ACKNOWLEDGEMENTS
%
% The acknowledgements environment puts a large, bold, centered 
% "Acknowledgements" label at the top of the page. The acknowledgements
% themselves appear in a quote environment, i.e. tabbed in at both sides, and
% on its own page.

\newenvironment{acknowledgements}
{\pagestyle{empty}
\addcontentsline{toc}{chapter}{Acknowledgements}
\begin{alwayssingle}
\begin{center}
\vspace*{1.5cm}
{\Large \bfseries Acknowledgements}
\end{center}
\vspace{0.5cm}
\begin{quote}}
{\end{quote}\end{alwayssingle}}

% The acknowledgementslong environment puts a large, bold, centered 
% "Acknowledgements" label at the top of the page. The acknowledgement itself 
% does not appears in a quote environment so you can get more in.

\newenvironment{acknowledgementslong}
{\pagestyle{empty}
\addcontentsline{toc}{chapter}{Acknowledgements}
\begin{alwayssingle}
\begin{center}
\vspace*{1.5cm}
{\Large \bfseries Acknowledgements}
\end{center}
\vspace{0.5cm}\begin{quote}}
{\end{quote}\end{alwayssingle}}

%ABSTRACT
%
%The abstract environment puts a large, bold, centered "Abstract" label at
%the top of the page. The abstract itself appears in a quote environment,
%i.e. tabbed in at both sides, and on its own page.

\newenvironment{abstracts} {\begin{alwayssingle} \pagestyle{empty}
\addcontentsline{toc}{chapter}{Abstract}
  \begin{center}
  \vspace*{1.5cm}
  {\Large \bfseries  Abstract}
  \end{center}
  \vspace{0.5cm}
   \begin{quote}}
{\end{quote}\end{alwayssingle}}

%The abstractlong environment puts a large, bold, centered "Abstract" label at
%the top of the page. The abstract itself does not appears in a quote
%environment so you can get more in.

\newenvironment{abstractslong} {\begin{alwayssingle} \pagestyle{empty}
\addcontentsline{toc}{chapter}{Abstract}
  \begin{center}
  \vspace*{1.5cm}
  {\Large \bfseries  Abstract}
  \end{center}
  \vspace{0.5cm} \begin{quote}}
{\end{quote}\end{alwayssingle}}

%The abstractseparate environment is for running of a page with the abstract
%on including title and author etc as required to be handed in separately

\newenvironment{abstractseparate} {\begin{alwayssingle} \pagestyle{empty}
  \vspace*{-1in}
 \begin{center}
    { \Large {\bfseries {\@title}} \par}
    {{\large \vspace*{1ex} \@author} \par}
{\large \vspace*{1ex}
    {{\@collegeordept} \par}
    {{\@university} \par}
\vspace*{1ex}
    {{\it \submittedtext} \par}
    {\it {\@degree} \par}
\vspace*{2ex}
    {\@degreedate}}
  \end{center}}
{\end{alwayssingle}}

%Statement of originality if required

\newenvironment{declaration} {\begin{alwayssingle} \pagestyle{empty}
\addcontentsline{toc}{chapter}{Declaration of Originality}
  \begin{center}
  \vspace*{1.5cm}
  {\Large \bfseries  Declaration of Originality}
  \end{center}
  \vspace{0.5cm}
   \begin{quote}}
{\end{quote}\end{alwayssingle}}


\newcommand{\copyrightpage}{
	\newpage
	\thispagestyle{empty}
	\vspace*{25pt}
	\begin{center}
	\scshape \noindent \small \copyright \  \small  \@author \\
	all rights reserved, \@degreeyear
	\end{center}
	\newpage
	\rm
}

%:-------------------------- page numbers: roman+arabic -----------------------

% ROMANPAGES
%
% The romanpages environment set the page numbering to lowercase roman one
% for the contents and figures lists. It also resets
% page-numbering for the remainder of the dissertation (arabic, starting at 1).

%\newenvironment{romanpages}
%{
%	\setcounter{page}{1}
%	\renewcommand{\thepage}{\roman{page}}
%} % close romanpage env't
 
{\newpage\renewcommand{\thepage}{\arabic{page}}\setcounter{page}{1}}





