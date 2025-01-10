# Check AMM

[![URL](https://img.shields.io/badge/URL-check--amm.vercel.app-blue)](https://check-amm.vercel.app/home)
[![DOI](https://img.shields.io/badge/DOI-10.1021%2Facs.orglett.4c04730-blue)](https://doi.org/10.1021/acs.orglett.4c04730)
![logo](https://github.com/user-attachments/assets/3edd8e9a-5b25-4905-b596-bde023825503)


## What is it?

`check-amm` is a Python-based software that processes Supporting Information (SI) documents to:

* Verify the integrity</strong> of Accurate Mass Measurement (AMM) data, formerly known as High-Resolution Mass Spectrometry (HRMS)
* Identify discrepancies</strong> between the calculated and measured accurate masses and the proposed molecular formula, and
* Comment on potential deviations</strong> from the proposed molecular formula that would reconcile the calculated and measured accurate masses.


## Installation and getting started

Coming soon...


## Usage

### Web App

The easiest way of utilizing this tool is through the [Check AMM website](https://check-amm.vercel.app/home). Instructions on how to operate the web app can be found on the website's [Help page](https://check-amm.vercel.app/help).

### Command Line

Coming soon...


## Support and Community

If you have questions, comments, or suggestions, the best place is [GitHub discussions](https://github.com/kozlowski-lab/check-amm/discussions). If you have found a bug or would like to request a feature, please [create an issue](https://github.com/kozlowski-lab/check-amm/issues).


## More about Check AMM

Check AMM arose from concerns about the lack of available tools for scrutinizing data quality in scientific publications. Research relies heavily on reproducibility; therefore, inconsistent or invalid data can hinder another researcher's ability to replicate findings, impeding scientific progress. Conclusions drawn from faulty data can be misleading and may propagate errors through subsequent research. These issues underscore the importance of ensuring high data quality, particularly when sharing work among the scientific community.

The Supporting Information (SI) of a scientific publication contains a combination of technical details and vast amounts of data. Due to the size of these documents, manual identification of all inaccuracies is nearly impossible for human reviewers. Moreover, there is a significant lack of tools for automated and standardized data quality assessments. To address this need, [Prof. Christmann](https://www.bcp.fu-berlin.de/en/chemie/chemie/forschung/OrgChem/christmann/index.html) from Freie Universität Berlin designed the [`HRMS-Checker`](https://github.com/match22lab/HRMS-Checker-2.0) program, which assesses the validity and internal consistency of Accurate Mass Measurement (AMM) data, formerly known as High-Resolution Mass Spectrometry (HRMS) in scientific publications. Inaccurate mass data can lead to incorrect compound identification or structural elucidation. In collaboration with the [Kozlowski lab](https://www.mckgroup.org/) from the University of Pennsylvania, this program, now named Check AMM, was adapted and made available to the community through this repository and the Check AMM web app.

The purpose of Check AMM is to provide scientists with a report assessing the validity and internal consistency of their accurate mass data before submitting their work for publication. This report offers an opportunity to identify and correct any inaccuracies and serves as a standard measure of AMM data reliability for reviewers and the scientific community.


## Acknowledgements

Thanks to the following contributors from the University of Pennsylvania:
- [Guillermo Correa Otero](https://github.com/guille797)
- [Sarah Zhang](https://github.com/sarahzhanng)

Special thanks to [Prof. Christmann](https://github.com/match22lab) for providing us with the original code to be adapted.


## Citation

How to cite:

Kozlowski, M. C.; Correa Otero, G.; Zhang, S. On the Integrity of Accurate Mass Measurement Data in Compound Characterization. *Org. Lett.* **2025**, *27* (1), 1–3. https://doi.org/10.1021/acs.orglett.4c04730.


## License

This project is licensed under the GNU GPLv3 License - see the [COPYING](COPYING) file for details. It also uses the following third-party libraries:
- [`molmass`](https://github.com/cgohlke/molmass) (BSD-3-Clause License)
- [`PyMuPDF`](https://github.com/pymupdf/PyMuPDF) (AGPL-3.0 License)
- [`pandas`](https://github.com/pandas-dev/pandas) (BSD-3-Clause license)
- [`FPDF`](https://github.com/Setasign/FPDF?tab=readme-ov-file) (Custom License)

The adapted [`HRMS-Checker`](https://github.com/match22lab/HRMS-Checker-2.0) is also licensed under the GNU GPLv3 License.
