# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a LaTeX thesis template for the University of Sydney (SAIL - Sydney Astrophotonics Instrumentation Laboratory). It's designed for PhD theses in astrophotonics and follows the university's thesis formatting requirements.

## Commands

### Compilation Commands
- `pdflatex thesis.tex` - Compile the main thesis document
- `biber thesis` - Process bibliography references using Biber
- `makeglossaries thesis` - Generate glossary and acronyms
- `pdflatex thesis.tex` (run twice after biber/makeglossaries for proper cross-references)

Full build sequence:
```bash
pdflatex thesis.tex
biber thesis
makeglossaries thesis
pdflatex thesis.tex
pdflatex thesis.tex
```

### Utility Scripts
- `sh tools/inkscapeSVG2PDFtex.sh <filename>` - Convert SVG to PDF with LaTeX text overlay
- `python tools/pdfcolorsplit <file.pdf>` - Split PDF into color and black/white sections for printing

## Architecture and Structure

### Document Structure
- `thesis.tex` - Main document that includes all chapters and sets up document structure
- `preamble.tex` - Contains all package imports, custom commands, title page setup, and styling
- `Latex/Classes/PhDthesisPSnPDF.cls` - Custom document class based on standard book class

### Content Organization
- `0_frontmatter/` - Abstract, acknowledgments, declarations, list of publications
- `0_overview/` - Template usage guide and documentation chapter
- `1_introduction/` - Introduction chapter
- `2/` through `7/` - Main thesis chapters (astrophotonics, design principles, spectrographs, etc.)
- `9_backmatter/` - Appendices, attached papers, glossary definitions
- `bib/` - Bibliography files (BibTeX format, managed with Biber)

### LaTeX Infrastructure
- `Latex/Classes/` - Custom document classes and bibliography styles
- `Latex/StyleFiles/` - Additional style packages (mcode.sty for MATLAB code, watermark.sty)
- `Latex/Macros/MacroFile1.tex` - Custom macros and convenience commands

### Key Custom Macros (from MacroFile1.tex)
- `\figuremacroW{filename}{title}{description}{width}` - Insert figure with specified width
- `\figuremacroSUB{title}{file1}{caption1}{width1}{file2}{caption2}{width2}` - Subfigure layout
- `\fig{label}`, `\Fig{label}` - Figure references with consistent formatting
- `\eqn{label}`, `\sect{label}`, `\chap{label}` - Cross-reference shortcuts
- `\includesvg{filename}` - Include SVG with automatic PDF conversion
- `\insertpaper{title}{label}{filename}` - Include full papers as PDF sections

### Draft vs Final Mode
The template supports draft mode controlled by:
- `\thesisisdrafttrue` - Enables draft mode with TODO lists and draft watermarks
- `\thesisisdraftfalse` - Final mode for submission

In draft mode, the template shows:
- List of TODOs on the first page
- Draft stamps in headers
- Comment system with color-coded authors (CHB, SLS, JGR)

### Bibliography System
- Uses Biber (not BibTeX) with biblatex
- References stored in `bib/thesisRefs.bib`
- Configured for segmented bibliography with `\bibbysegment`
- Bibliography style optimized for academic physics/astronomy citations

### Graphics and Figures
- Graphics paths automatically set for all chapter figure directories
- Supports automatic SVG to PDF conversion with LaTeX text overlay
- Inkscape integration for vector graphics workflow
- Figure macros provide consistent formatting and labeling

### Special Features
- **Template Overview Chapter** (`0_overview/overview.tex`) - Comprehensive guide explaining how to use the template, customize metadata, manage bibliography, use figure macros, and implement best practices
- Glossary and acronym management with `makeglossaries`
- Code listings optimized for MATLAB with syntax highlighting
- University branding with custom colors and logos
- Comment system for collaborative editing with different author colors
- Tools for print optimization (color/B&W splitting)
- Draft mode support for faster compilation during development

### Chapter Management
Use `\includeonly{}` commands in thesis.tex to compile specific chapters during development. Examples are commented out in the main file for each chapter.

## Important Notes

- This template is specifically configured for University of Sydney thesis requirements
- Uses custom PhD thesis class that extends the book class
- Bibliography must be processed with Biber, not BibTeX
- Graphics workflow assumes Inkscape for SVG processing
- Final compilation requires multiple LaTeX runs for proper cross-references