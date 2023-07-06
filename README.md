# LeadCrystal
Code and data from measurement of Lead (Pb) crystal leaching lead into beverages

Link to Google Docs spreadsheet (same as the xlsx in this repo)
https://docs.google.com/spreadsheets/d/1cfIEGUUV5k4__KO8_MILYz-mCUhh7hRl3U3gHKuC_7I/edit?usp=sharing

Link to video showcasing the whole project
https://youtu.be/-cAB5FG4bXI

This project makes use of dithizone, a chemical that produces colorful products when mixed with trace amounts of metals. https://ca.hach.com/asset-get.download.jsa?id=7639983726
5ml dichloromethane was used to extract the lead from 30ml aqueous samples containing 0.5% acetic acid, after adjusting pH with NH4OH to 8.7.  The color of the sample was analyzed in 10mm glass cuvettes with a LR1-B spectrometer
https://www.aseq-instruments.com/LR1.html
Exposure was 15ms, illumination provided by an 18 watt tungsten lamp

The main DCM extraction method
https://ca.hach.com/asset-get.download.jsa?id=7639983726
https://archive.org/details/gov.law.apha.method.3500-pb.1992/page/n5/mode/2up
https://pubs.usgs.gov/bul/1084f/report.pdf



https://sci-hub.se/https://doi.org/10.2116/analsci.23.193
Us4876068a

https://www.lenntech.com/periodic/water/lead/lead-and-water.htm


Lead cystal study.xslx - summary spreadsheet of results

Pb_calibration.m - Octave script to collect data, and generate calibration curve based on standards of lead

Pb_sample.m - Octave script to collect data, and use calibration curve to determine concentrations

calfit - Octave data file containing only the required calibration coeficients to determine lead concentration in unknown samples

calsession - the whole Octave workspace used during calibration determination
samplesession - the whole Octave workspace used when calculating unknown lead concentrations
caldata - Octave variables used during development -- not needed for current operation
