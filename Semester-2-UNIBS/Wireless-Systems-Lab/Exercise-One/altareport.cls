%%%%%%%%%%%%%%%%%
% This is altareport.cls (v0.0.1, 24 April 2021) written by
% Andrew Wilson (andrew.s.wilson@tutanota.com).
%
%% It may be distributed and/or modified under the
%% conditions of the LaTeX Project Public License, either version 1.3
%% of this license or (at your option) any later version.
%% The latest version of this license is in
%%    http://www.latex-project.org/lppl.txt
%% and version 1.3 or later is part of all distributions of LaTeX
%% version 2003/12/01 or later.
%%
%%
% Contributions:
% - None... yet.

%%%%%%%%%%%%%%%%
\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{altareport}[2021/04/24 - Adapted from AltaCV class for resumes.]

%% v1.1: Optionally load academicons
\newif\if@academicons
\DeclareOption{academicons}{\@academiconstrue}
%% v1.1.3: Choice of round/square photo
\newif\if@normalphoto
\DeclareOption{normalphoto}{\@normalphototrue}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{extarticle}}
\newif\if@raggedtwoe
\DeclareOption{ragged2e}{\@raggedtwoetrue}
%% v1.3: load hyperref for clickable hyperlinks
\newif\if@withhyper
\DeclareOption{withhyper}{\@withhypertrue}
\ProcessOptions\relax

\LoadClass{article}
%% v1.3.2 Hopefully this helps make the PDF
%% file more 'friendly' with copy-paste etc
\RequirePackage[a-1b]{pdfx}
%% v1.3.3 Temporary fix for pdfx and biblatex incompatibility
%% (Thanks to @badetitou for highlighting this issue)
%% https://github.com/plk/biblatex/issues/1044
%% Will be removed when fix is released in biblatex 3.16
\catcode30=12
\RequirePackage{geometry}
\RequirePackage[fixed]{fontawesome5}
\RequirePackage{scrlfile}
\RequirePackage{xparse}



%% v1.1: Optionally load academicons
%% v1.1.5: Handle different versions of academicons
\if@academicons
    \ClassError{altareport}{academicons unsupported by latex or pdflatex.}
\fi

\if@raggedtwoe
  \RequirePackage[document]{ragged2e}
\fi

\if@withhyper
  \AtBeginDocument{%
    \RequirePackage{hyperref}
    \hypersetup{
    	linkcolor=body,
    	citecolor=accent,
    	linktoc=all,
    	urlcolor=accent}

    \urlstyle{
    	colorlinks=true,
    	citecolor=accent,
  		linkcolor=body,
  		linktoc=all,
  		urlcolor=accent}
  	}
\fi


\RequirePackage{xcolor} %% Pretty Colours
\RequirePackage{tikz}
\usetikzlibrary{arrows}
\RequirePackage[skins]{tcolorbox}
\RequirePackage[inline]{enumitem}
\setlist{leftmargin=*,labelsep=0.5em,nosep,itemsep=0.25\baselineskip,after=\vspace{0.25\baselineskip}}
\setlist[itemize]{label=\itemmarker}
\RequirePackage{graphicx} %% Pictures, PDFs, etc.
\RequirePackage{etoolbox}
\RequirePackage{dashrule}
\RequirePackage{multirow,tabularx}
\RequirePackage{changepage}
% \RequirePackage{marginfix}
\RequirePackage{colortbl}
\RequirePackage{caption}

\newcommand{\itemmarker}{{\small\textbullet}}
\newcommand{\ratingmarker}{\faCircle}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% DATE FORMAT %%%%
\RequirePackage[nodayofweek]{datetime} % Auto update the date when file is compiled
\newdateformat{monthyeardate}{\monthname[\THEMONTH], \THEYEAR} % This will produce the month, year format when \monthname\today is used.


%%%%%%%%%% THEMES %%%%%%%%%%

\newcommand{\PastelBlue}{
	\definecolor{SlateGrey}{HTML}{2E2E2E}
	\definecolor{LightGrey}{HTML}{4C4C4C}
	\definecolor{DarkerPastelTurquoise}{HTML}{00535B}
	\definecolor{DarkPastelTurquoise}{HTML}{006B74}
	\definecolor{PastelTurquoise}{HTML}{009AA7}
	\definecolor{Silver}{HTML}{C0C0C0}

	\colorlet{heading}{DarkerPastelTurquoise}
	\colorlet{headingrule}{Silver}
	\colorlet{subheading}{DarkPastelTurquoise}
	\colorlet{name}{DarkPastelTurquoise}
	\colorlet{accent}{PastelTurquoise}
	\colorlet{tagline}{PastelTurquoise}
	\colorlet{emphasis}{SlateGrey}
	\colorlet{body}{LightGrey}
}


\newcommand{\GreenAndGold}{
	\definecolor{SlateGrey}{HTML}{2E2E2E}
	\definecolor{LightGrey}{HTML}{4C4C4C}
	\definecolor{PastelGreen}{HTML}{00a86b}
	\definecolor{DarkerPastelGreen}{HTML}{00422a}
	\definecolor{DarkPastelGreen}{HTML}{00754b}
	\definecolor{DarkGoldenRod}{HTML}{b8860b}

	\colorlet{heading}{DarkerPastelGreen}
	\colorlet{headingrule}{DarkGoldenRod}
	\colorlet{subheading}{DarkPastelGreen}
	\colorlet{name}{DarkPastelGreen}
	\colorlet{accent}{PastelGreen}
	\colorlet{tagline}{PastelGreen}
	\colorlet{emphasis}{SlateGrey}
	\colorlet{body}{LightGrey}
}

\newcommand{\Purple}{
	\definecolor{SlateGrey}{HTML}{2E2E2E}
	\definecolor{LightGrey}{HTML}{4C4C4C}
	\definecolor{DarkerMagenta}{HTML}{4C344F}
	\definecolor{DarkMagenta}{HTML}{69486E}
	\definecolor{Magenta}{HTML}{875C8D}
	\definecolor{GoldenRod}{HTML}{DAA520}
	\definecolor{Silver}{HTML}{C0C0C0}

	\colorlet{heading}{DarkerMagenta}
	\colorlet{headingrule}{GoldenRod}
	\colorlet{subheading}{DarkMagenta}
	\colorlet{name}{DarkMagenta}
	\colorlet{accent}{Magenta}
	\colorlet{tagline}{Magenta}
	\colorlet{emphasis}{SlateGrey}
	\colorlet{body}{LightGrey}
}

\newcommand{\UNILIM}{
	\definecolor{SlateGrey}{HTML}{2E2E2E}
	\definecolor{LightGrey}{HTML}{4C4C4C}
	\definecolor{DarkerPastelRed}{HTML}{710E15}
	\definecolor{DarkPastelRed}{HTML}{881019}
	\definecolor{PastelRed}{HTML}{B51621}
	\definecolor{GoldenEarth}{HTML}{E7D192}
	\definecolor{Silver}{HTML}{C0C0C0}

	\colorlet{heading}{DarkerPastelRed}
	\colorlet{headingrule}{Silver}
	\colorlet{subheading}{DarkPastelRed}
	\colorlet{name}{DarkPastelRed}
	\colorlet{accent}{PastelRed}
	\colorlet{tagline}{PastelRed}
	\colorlet{emphasis}{SlateGrey}
	\colorlet{body}{LightGrey}
}

\newcommand{\UNIBS}{
	\definecolor{SlateGrey}{HTML}{2E2E2E}
	\definecolor{LightGrey}{HTML}{4C4C4C}
	\definecolor{DarkerModerateBlue}{HTML}{324879}
	\definecolor{DarkModerateBlue}{HTML}{3D5895}
	\definecolor{ModerateBlue}{HTML}{738ACC}
	\definecolor{SlightlyDesaturatedBlue}{HTML}{607DBF}
	\definecolor{Silver}{HTML}{A7AABD}

	\colorlet{heading}{DarkerModerateBlue}
	\colorlet{headingrule}{Silver}
	\colorlet{subheading}{DarkModerateBlue}
	\colorlet{name}{DarkModerateBlue}
	\colorlet{accent}{ModerateBlue}
	\colorlet{tagline}{ModerateBlue}
	\colorlet{emphasis}{SlateGrey}
	\colorlet{body}{LightGrey}
}

\newcommand{\PastelRed}{
	\definecolor{SlateGrey}{HTML}{2E2E2E}
	\definecolor{LightGrey}{HTML}{666666}
	\definecolor{DarkerPastelRed}{HTML}{600909}
	\definecolor{DarkPastelRed}{HTML}{8F0D0D}
	\definecolor{PastelRed}{HTML}{BD1111}
	\definecolor{GoldenEarth}{HTML}{E7D192}

	\colorlet{heading}{DarkerPastelRed}
	\colorlet{headingrule}{GoldenEarth}
	\colorlet{subheading}{DarkPastelRed}
	\colorlet{name}{DarkPastelRed}
	\colorlet{accent}{PastelRed}
	\colorlet{tagline}{PastelRed}
	\colorlet{emphasis}{SlateGrey}
	\colorlet{body}{LightGrey}
}

\newcommand{\BoringDefault}{
	\definecolor{SlateGrey}{HTML}{2E2E2E}
	\definecolor{LightGrey}{HTML}{4C4C4C}
	\definecolor{Black}{HTML}{000000}
	\definecolor{Grey}{HTML}{a9a9a9}

	\colorlet{heading}{Black}
	\colorlet{headingrule}{Grey}
	\colorlet{subheading}{SlateGrey}
	\colorlet{name}{SlateGrey}
	\colorlet{accent}{SlateGrey}
	\colorlet{tagline}{SlateGrey}
	\colorlet{emphasis}{SlateGrey}
	\colorlet{body}{LightGrey}
}


%%%%% ADD NEW THEMES BELOW %%%%%
\newcommand\SetTheme[1]{
	\ifstrequal{#1}{UNIBS}{\UNIBS}{ % else
		\ifstrequal{#1}{UNILIM}{\UNILIM}{ % else
			\ifstrequal{#1}{PastelBlue}{\PastelBlue}{ % else
				\ifstrequal{#1}{GreenAndGold}{\GreenAndGold}{ % else
  					\ifstrequal{#1}{PastelRed}{\PastelRed}{ % else
  						\ifstrequal{#1}{Purple}{\Purple}{\BoringDefault}
					}
				}
			}
		}
	}
}

%%%%%%%%%% HEADER & FOOTER %%%%%%%%%%
\RequirePackage{fancyhdr}

\fancyhead[L]{\color{subheading}\FootNoteFont{\leftmark}}
\fancyhead[C]{}
\fancyhead[R]{\color{LightGrey}\FootNoteFont{\@title}}
\fancyfoot[L]{\color{subheading}\FootNoteFont{\@author}}
\fancyfoot[C]{}
\fancyfoot[R]{\color{LightGrey}\FootNoteFont{\thepage}}
\renewcommand{\headrulewidth}{2pt}
\renewcommand{\footrulewidth}{0pt}

\renewcommand{\headrule}{\hbox to\headwidth{%
  \color{headingrule}\leaders\hrule height \headrulewidth\hfill}}


\setlength{\headheight}{20pt}% ...at least 18.6....pt
%%\setlength{\footheight}{20pt}% ...at least 18.6....pt


%%%%%%%%%% SECTIONS & SUBSECTIONS & SUBSUBSECTIONS %%%%%%%%%%
\RequirePackage{titlesec}

\titleformat{\section}[hang]
{\color{heading}\ReportSectionFont}
{\color{heading}\thesection}{5mm}{}{}

\titleformat{\subsection}[hang]
{\color{subheading}\ReportSubSectionFont}
{\color{subheading}\thesubsection}{5mm}{}{}

\titleformat{\subsubsection}[hang]
{\color{accent}\ReportEmphasisFont}
{\color{accent}\thesubsubsection}{5mm}{}{}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\setlength{\parindent}{0pt}
\newcommand{\divider}{\textcolor{body!30}{\hdashrule{\linewidth}{0.6pt}{0.5ex}}\medskip}

\newenvironment{fullwidth}{%
  \begin{adjustwidth}{}{\dimexpr-\marginparwidth-\marginparsep\relax}}
  {\end{adjustwidth}}

%% v1.3.1 \detokenize will break UTF-8 in pdflatex
%% Using alternative from https://tex.stackexchange.com/a/530911/226
\newcommand{\utffriendlydetokenize}[1]{%
\scantokens{%
  \catcode`\_=12%
%   \catcode`\^=12%
%   \catcode`\{=12%
%   \catcode`\}=12%
  \catcode`\&=12%
  \catcode`\$=12%
  \catcode`\#=12%
  \catcode`\~=12%
%   \catcode`\\=12%
  {#1}%
}%
}
%% v1.3: Incorporating hyperlinks
%% v1.3.1: using \unfriendlydetokenize to avoid
%% breaking unicode
\ExplSyntaxOn
\NewDocumentCommand{\printinfo}{m m o}{%
  \IfNoValueTF{#3}{%
    \mbox{\textcolor{accent}{\normalfont #1}~\utffriendlydetokenize{#2}\hspace{2em}}%
  }{%
     \if@withhyper%
      \mbox{\textcolor{accent}{\normalfont #1}~
      \href{#3}{\utffriendlydetokenize{#2}}\hspace{2em}}
     \else%
       \ClassWarning{Please specify [withhyper] option to enable hyperlinks. Printing out full hyperlink prefix #1 for now.}%
      \mbox{\textcolor{accent}{\normalfont #1}~{\utffriendlydetokenize{#3#2}}\hspace{2em}}%
     \fi%
  }%
}%

%% v1.3: Exploring convenient creation of fields
\NewDocumentCommand{\NewInfoField}{m m o}{%
  \IfNoValueF{#3}{\csdef{#1 hyperprefix}{#3}}%
  \csdef{#1 symbol}{#2}%
  \csdef{#1}##1{%
    \if@withhyper
      \IfNoValueTF {#3}
        {\printinfo{\csuse{#1 symbol}}{##1}}%
        {\printinfo{\csuse{#1 symbol}}{##1}[\csuse{#1 hyperprefix}##1]}%
    \else
      \printinfo{\csuse{#1 symbol}}{##1}%
    \fi%
  }
}
\ExplSyntaxOff


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% FONT DEFINITIONS %%%%%%%%%%
\newcommand{\ReportTitleFont}{\Huge\bfseries}
\newcommand{\ReportSubTitleFont}{\huge\bfseries}
\newcommand{\ReportSectionFont}{\LARGE\bfseries}
\newcommand{\ReportSubSectionFont}{\large\bfseries}
\newcommand{\ReportEmphasisFont}{\bfseries}
\newcommand{\FootNoteFont}{\footnotesize\bfseries}

\DeclareCaptionFont{CaptionLabel}{\color{body}\ReportEmphasisFont}
\DeclareCaptionFont{CaptionText}{\color{body}}
\captionsetup{labelfont={CaptionLabel},textfont={CaptionText}}


%%%%%%%%%% TITLE PAGE %%%%%%%%%%
%%%%%%%%%% TITLE PAGE VARIABLES %%%%%%%%%%
\newcommand{\ReportTitle}[1]{\def\@title{#1}}
\newcommand{\SubTitle}[1]{\def\@subtitle{#1}}
\newcommand{\Author}[1]{\def\@author{#1}}
\newcommand{\ModCoord}[1]{\def\@coord{#1}}
\newcommand{\ReportDate}[1]{\def\@reportdate{#1}}
\newcommand{\FacultyOrLocation}[1]{\def\@location{#1}}


\newcommand{\MakeReportTitlePage}{
	\thispagestyle{empty}
	\vfill
	\vspace*{\fill}
	\begin{center}
	
	{\color{heading}\ReportTitleFont{\@title}}\bigskip \\
	{\color{subheading}\ReportSubTitleFont{\@subtitle}}\medskip \\
	{\bfseries{\@author}\smallskip} \\
	\@reportdate
	
	\vspace{2cm}
	\includegraphics[width=10cm]{Images/University-of-Brescia-Logo.png}
	\vspace{2cm}

	\large{\@location \\
	Module Coordinator: \@coord} \\

	\end{center}
	\vspace*{\fill}
	\vfill
	\newpage
	\pagestyle{fancy}
}


%%%%%%%%%% PERSONAL INFO %%%%%%%%%%
%%%%%%%%%% PERSONAL INFO VARIABLES %%%%%%%%%%
\newcommand{\name}[1]{\def\@name{#1}}
\newcommand{\tagline}[1]{\def\@tagline{#1}}
\newcommand{\personalinfo}[1]{\def\@personalinfo{#1}}

\NewInfoField{email}{\faAt}[mailto:]
\NewInfoField{mailaddress}{\faEnvelope}
\NewInfoField{phone}{\faPhone}[tel:]
\NewInfoField{homepage}{\faGlobe}[https://]
\NewInfoField{twitter}{\faTwitter}[https://twitter.com/]
\NewInfoField{linkedin}{\faLinkedin}[https://linkedin.com/in/]
\NewInfoField{github}{\faGithub}[https://github.com/]
\NewInfoField{orcid}{\aiOrcid}[https://orcid.org/]
\NewInfoField{location}{\faMapMarker}


\newcommand{\namefont}{\Large\bfseries}
\newcommand{\taglinefont}{\normalsize\bfseries}
\newcommand{\personalinfofont}{\footnotesize\bfseries}


\newcommand{\makeauthordetails}{%
	\begingroup
		\altacv@left@photos\hfill%
		\begin{minipage}{\dimexpr\linewidth-\altacv@photos@width}%
			\raggedright%
			{\namefont\color{name}{\@name}\par}
			\medskip
			{\taglinefont\color{tagline}\@tagline\par}
			\smallskip
			{\personalinfofont\@personalinfo\par}
			\bigskip
		\end{minipage}\hfill%
		\altacv@right@photos\par%
		{\color{headingrule}\rule{\linewidth}{2pt}\par}
	\endgroup\medskip
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% v1.2: Support for multiple photos
\newlength{\altacv@photos@width}
\newlength{\altacv@photo@diam@left}
\newlength{\altacv@photo@diam@right}
\def\altacv@left@photos{}
\def\altacv@right@photos{}

\newcommand{\@makeaphoto}[2]{%
  \begin{minipage}{#1}%
    \if@normalphoto
      \includegraphics[width=\linewidth]{#2}
    \else
      \tikz\path[fill overzoom image={#2}]circle[radius=0.5\linewidth];
    \fi%
  \end{minipage}%
}

\newcommand{\altacv@add@photo@left}[1]{%
  \appto{\altacv@left@photos}{%
    \@makeaphoto{\altacv@photo@diam@left}{#1}\hspace{1ex}%
  }%
  \addtolength{\altacv@photos@width}{\altacv@photo@diam@left}%
  \addtolength{\altacv@photos@width}{1ex}%
}
\newcommand{\altacv@add@photo@right}[1]{%
  \appto{\altacv@right@photos}{%
    \@makeaphoto{\altacv@photo@diam@right}{#1}\hspace{1ex}%
  }%
  \addtolength{\altacv@photos@width}{\altacv@photo@diam@right}%
  \addtolength{\altacv@photos@width}{1ex}%
}
\newcommand{\photoL}[2]{%
  \setlength{\altacv@photo@diam@left}{#1}%
  \forcsvlist{\altacv@add@photo@left}{#2}%
}
\newcommand{\photoR}[2]{%
  \setlength{\altacv@photo@diam@right}{#1}%
  \forcsvlist{\altacv@add@photo@right}{#2}%
}
\let\photo\photoR




\renewenvironment{quote}{\color{accent}\itshape\large}{\par}

\newenvironment{reportcolumn}[1]{\begin{minipage}[t]{#1}\raggedright}{\end{minipage}}


% v1.1.2: make it easier to add a sidebar aligned with top of next page
\RequirePackage{afterpage}
\newcommand{\addsidebar}[2][]{\marginpar{%
  \ifstrequal{#1}{}{}{\vspace*{#1}}%
  \input{#2}}%
}
\newcommand{\addnextpagesidebar}[2][]{\afterpage{\addsidebar[#1]{#2}}}

\AtBeginEnvironment{tabular}{\color{body}\normalsize}

\AtBeginDocument{%
	\pagestyle{empty}
	\color{body}
	\raggedright
	\setlength\arrayrulewidth{0.5pt}
	\arrayrulecolor{body}
}
