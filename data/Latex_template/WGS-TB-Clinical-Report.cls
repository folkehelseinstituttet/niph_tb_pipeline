%% Copyright 2017 Anamaria Crisan <acrisan@cs.ub.ca>

%-------------------------------------
% Identification
%-------------------------------------
\ProvidesClass{WGS-TB-Clinical-Report}[2017/07/12 v1.0 WGS TB Clinical Report Design]
\NeedsTeXFormat{LaTeX2e}

%-------------------------------------
% Class options
%-------------------------------------

% Options for draft or final
\DeclareOption{draft}{\setlength\overfullrule{5pt}}
\DeclareOption{final}{\setlength\overfullrule{0pt}}
% Inherit options of article
\DeclareOption*{
  \PassOptionsToClass{\CurrentOption}{article}
}
\ProcessOptions\relax
\LoadClass{article}

%-------------------------------------
% 3rd party packages
%-------------------------------------
\usepackage{fancyhdr}
\usepackage{lastpage}

\newcommand*{\headerstyle}[1]{{\fontsize{18pt}}}

%-------------------------------------
% New Commands
%-------------------------------------

% Define the elements of the report header
% Usage: \makeReportHeader{<left>}{<right>}
\newcommand*{\makeReportHeader}[2]{
  \fancyhead{}
  \fancyhead[L]{
    \headerstyle{#1}
  }
  \fancyhead[C]{
    \headerstyle{#2}
  } 
}
