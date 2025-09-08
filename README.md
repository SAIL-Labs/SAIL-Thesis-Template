# SAIL Thesis Template

A comprehensive LaTeX thesis template based on Jakob Suckale's PhDthesisPSnPDF class, customized for the Uni Sydney requirmetns for Physics students and academic thesis writing.

## Features

- **Professional Layout**: Book-style document with proper margins and typography
- **Modular Structure**: Organized chapters and sections for easy maintenance
- **Bibliography Support**: Integrated biblatex with custom styles
- **Figure Management**: Dedicated figure directories with subfigure support macros.
- **Glossary Support**: Built-in glossary and acronym management, prefilled with some common stuff for astrophotonics.
- **Custom Macros**: some Astronomy-specific macros and shortcuts
- **PDF Optimization**: Hyperlinked table of contents and references

## Quick Start

1. **Clone or download** this template
2. **Edit the main content files** in each chapter directory
3. **Update the bibliography** in `bib/thesisRefs.bib`
4. **Customize frontmatter** (abstract, acknowledgments, etc.)
5. **Compile** using your preferred LaTeX editor (recommended: Overleaf)

## Directory Structure

```
SAIL-Thesis-Template/
├── README.md                   # This documentation
├── thesis.tex                  # Main document file
├── preamble.tex               # Document preamble
├── CLAUDE.md                  # Claude Code configuration
│
├── 0_frontmatter/             # Front matter components
│   ├── abstract.tex           # Thesis abstract
│   ├── acknowledgement.tex    # Acknowledgments
│   ├── declaration.tex        # Declaration of originality
│   ├── dedication.tex         # Dedication page
│   ├── list_of_publications.tex # Publication list
│   └── figures/               # Frontmatter figures
│
├── 1_introduction/            # Introduction chapter
│   ├── introduction.tex       # Main introduction content
│   ├── overview.tex          # Overview section
│   └── figures/              # Introduction figures
│
├── 2/                        # Chapter 2 directory
│   ├── macros.tex            # Chapter-specific macros
│   └── figures/              # Chapter 2 figures
│
├── 3/                        # Chapter 3 directory
│   ├── figures.tex           # Chapter content
│   └── figures/              # Chapter 3 figures
│
├── 9_backmatter/             # Back matter components
│   ├── appendix.tex          # Appendices
│   ├── glossary.tex          # Glossary definitions
│   ├── attached_papers.tex   # Attached papers section
│   ├── python_example.py     # Python code example
│   ├── matlab_example.m      # MATLAB code example
│   ├── zemax_example.zpl     # ZEMAX ZPL macro example
│   └── papers/               # Attached paper files
│
├── bib/                      # Bibliography
│   └── thesisRefs.bib        # BibTeX references
│
├── Latex/                    # LaTeX configuration
│   ├── Classes/              # Document classes
│   │   └── PhDthesisPSnPDF.cls # Main thesis class
│   ├── Macros/               # Custom macros
│   │   ├── MacroFile1.tex    # General macros
│   │   └── aas_macros.sty    # Astronomy macros
│   └── StyleFiles/           # Style files
│       ├── mcode.sty         # MATLAB code formatting
│       └── watermark.sty     # Watermark support
│
└── tools/                    # Additional tools and utilities
```

## File Organization

### Naming Convention
- **Frontmatter**: Prefix `0_` (e.g., `0_frontmatter/`)
- **Main chapters**: Numbered directories (e.g., `1_introduction/`, `2/`, `3/`)
- **Backmatter**: Prefix `9_` (e.g., `9_backmatter/`)

### Chapter Structure
Each chapter directory should contain:
- Main content file (`.tex`)
- `figures/` subdirectory for images
- Optional chapter-specific macros or includes

## Frontmatter Components

### Abstract (`0_frontmatter/abstract.tex`)
- Uses `\begin{abstracts}...\end{abstracts}` environment
- Alternative `\begin{abstractslong}...\end{abstractslong}` for longer abstracts
- Currently uses placeholder text (`\lipsum`)

### Declaration (`0_frontmatter/declaration.tex`)
- Statement of originality and authorship
- Includes space for signature (electronic signature supported)
- Customizable for different university requirements

### Other Frontmatter
- **Acknowledgments**: Thank supervisors, colleagues, funding sources
- **Dedication**: Optional personal dedication
- **List of Publications**: Your published works related to the thesis

## Main Content Chapters

### Chapter 1: Introduction (`1_introduction/`)
- `introduction.tex`: Main introduction content
- `overview.tex`: Thesis overview and structure
- Establishes context and research questions

### Subsequent Chapters (`2/`, `3/`, etc.)
- Add new numbered directories for each chapter
- Include main `.tex` file and `figures/` directory
- Use consistent naming and structure

## Backmatter Components

### Appendices (`9_backmatter/appendix.tex`)
Currently includes:

#### Code Examples
- **Python Class Example**: Optical system analysis class
- **MATLAB Class Example**: Wavefront analyzer class  
- **ZEMAX ZPL Macro**: Optical system automation script

#### Custom Language Support
- ZEMAX ZPL syntax highlighting defined with `\lstdefinelanguage{zpl}`
- Extensive keyword list for proper code formatting

### Glossary (`9_backmatter/glossary.tex`)
- Define technical terms and acronyms
- Uses `glossaries` package for automatic management
- Generates sorted glossary and acronym lists

### Bibliography
- Located in `bib/thesisRefs.bib`
- Uses biblatex with authoryear-comp style
- Automatically generates "References" section

## LaTeX Configuration

### Document Class (`Latex/Classes/PhDthesisPSnPDF.cls`)
Based on Jakob Suckale's PhD thesis class with features:
- Book-style layout with proper margins
- Custom chapter and section formatting
- Integrated hyperlinks and PDF bookmarks
- Bibliography and glossary support

### Key Packages Included
- `amsmath`, `amssymb`: Mathematical symbols and equations
- `graphicx`: Figure inclusion and manipulation
- `biblatex`: Advanced bibliography management
- `glossaries`: Glossary and acronym support
- `listings`: Source code formatting
- `hyperref`: PDF hyperlinks and bookmarks
- `caption`, `subcaption`: Figure caption control

### Custom Macros
- **General macros**: `Latex/Macros/MacroFile1.tex`
- **Astronomy macros**: `Latex/Macros/aas_macros.sty`
- Chapter-specific macros in individual directories

## Code Listings

### Supported Languages
- **Python**: Full syntax highlighting
- **MATLAB**: Using `mcode.sty` package
- **ZEMAX ZPL**: Custom language definition with extensive keywords

### Usage Example
```latex
\lstinputlisting[caption={Example Python code},
                 language=Python,
                 basicstyle=\tiny,
                 label={code:python_example}]{9_backmatter/python_example.py}
```

## Building the Document

### Recommended: Overleaf
1. Upload template to Overleaf
2. Set main document to `thesis.tex`
3. Compile using XeLaTeX or pdfLaTeX
4. Bibliography will be processed automatically

### Local Compilation
1. Ensure full LaTeX distribution (TeXLive/MikTeX)
2. Compile sequence:
   ```bash
   pdflatex thesis
   biber thesis
   makeglossaries thesis
   pdflatex thesis
   pdflatex thesis
   ```

### Required Packages
Most packages are standard, but ensure you have:
- `biblatex` with biber backend
- `glossaries` with makeglossaries
- `listings` for code highlighting
- Modern LaTeX distribution (2018+)

## Customization

### Adding New Chapters
1. Create new numbered directory (e.g., `4/`)
2. Add main content file (e.g., `chapter4.tex`)
3. Create `figures/` subdirectory
4. Include chapter in main `thesis.tex` file

### Modifying Layout
- Edit `Latex/Classes/PhDthesisPSnPDF.cls` for document class changes
- Adjust margins, fonts, and spacing in class file
- Custom headers/footers via `fancyhdr` package

### Bibliography Style
- Change style in class file: `\RequirePackage[style=authoryear-comp]{biblatex}`
- Supported styles: `numeric`, `alphabetic`, `authoryear`, `authortitle`
- Custom .bst files can be added to `Latex/StyleFiles/`

### Adding New Code Languages
1. Define language in appendix.tex:
```latex
\lstdefinelanguage{mylang}{
    morekeywords={keyword1, keyword2},
    morecomment=[l]{//},
    morestring=[b]"
}
```
2. Use with `\lstinputlisting[language=mylang]{file.ext}`

## Troubleshooting

### Common Issues

#### "File not found" errors
- Check file paths are correct and relative to main document
- Ensure figure files exist in specified directories
- Verify bibliography file path in biblatex setup

#### Bibliography not appearing
- Run biber after initial compilation
- Check .bib file syntax for errors
- Ensure citations exist in document (`\cite{}`)

#### Glossary not generating
- Run makeglossaries between compilations
- Check glossary entries are properly defined
- Ensure glossary package is loaded correctly

#### Code listings not formatting
- Verify language is supported or properly defined
- Check file encoding (UTF-8 recommended)
- Ensure listings package is loaded

### Best Practices

#### Version Control
- Use git to track changes
- Commit frequently with descriptive messages
- Consider using .gitignore for auxiliary LaTeX files

#### Figure Management
- Use vector formats (PDF, EPS) when possible
- Keep original high-resolution figures separate
- Use consistent naming convention
- Include alt text in captions for accessibility

#### Writing Workflow
1. Start with outline in comments
2. Write content in small sections
3. Compile frequently to catch errors
4. Use consistent citation style throughout
5. Proofread final version carefully

## Contributing

### Reporting Issues
- Check existing documentation first
- Provide minimal working example of problems
- Include LaTeX log output for compilation errors

### Improvements
- Fork the repository
- Make changes in feature branch
- Test thoroughly before submitting
- Document any new features or changes

### Template Updates
- Maintain backward compatibility when possible
- Update documentation for any changes
- Test with multiple LaTeX distributions
- Consider accessibility and internationalization

## License and Acknowledgments

This template is based on:
- **PhDthesisPSnPDF v2** by Jakob Suckale (2007)
- **CUEDthesis v1** by Harish Bhanderi (2002)

Customized for Sydney Institute for Astronomy (SAIL) use cases with additional features for scientific thesis writing.

---

For questions, issues, or suggestions, please refer to the documentation or contact the template maintainers.