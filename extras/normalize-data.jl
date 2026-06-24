# # Normalize CaMKAR data
using CSV
using DataFrames

datadir = joinpath(pwd(), "docs/data")
# ## Frequency data
freqdf = CSV.read(joinpath(datadir, "CaMKAR-freq.csv"), DataFrame)

onehz_initial = freqdf[!, "1Hz (Mean)"][1]
twohz_initial = freqdf[!, "2Hz (Mean)"][1]

# Normalize the mean and SD relative to the initial value
freqdf[!, "1Hz (Mean)"] .= freqdf[!, "1Hz (Mean)"] ./ onehz_initial
freqdf[!, "1Hz (SD)"] .= freqdf[!, "1Hz (SD)"] ./ onehz_initial
freqdf[!, "2Hz (Mean)"] .= freqdf[!, "2Hz (Mean)"] ./ twohz_initial
freqdf[!, "2Hz (SD)"] .= freqdf[!, "2Hz (SD)"] ./ twohz_initial

# Save the normalized data to a new CSV file
CSV.write(joinpath(datadir, "CaMKAR-freq-normalized.csv"), freqdf)

# ## Duration data
durationdf = CSV.read(joinpath(datadir, "CaMKAR-duration.csv"), DataFrame)
fifteen_initial = durationdf[!, "1Hz 15sec (Mean)"][1]
thirty_initial = durationdf[!, "1Hz 30sec (Mean)"][1]
sixty_initial = durationdf[!, "1Hz 60sec (Mean)"][1]
ninety_initial = durationdf[!, "1Hz 90sec (Mean)"][1]

# Normalize the mean and SD relative to the initial value
durationdf[!, "1Hz 15sec (Mean)"] .= durationdf[!, "1Hz 15sec (Mean)"] ./ fifteen_initial
durationdf[!, "1Hz 15sec (SD)"] .= durationdf[!, "1Hz 15sec (SD)"] ./ fifteen_initial
durationdf[!, "1Hz 30sec (Mean)"] .= durationdf[!, "1Hz 30sec (Mean)"] ./ thirty_initial
durationdf[!, "1Hz 30sec (SD)"] .= durationdf[!, "1Hz 30sec (SD)"] ./ thirty_initial
durationdf[!, "1Hz 60sec (Mean)"] .= durationdf[!, "1Hz 60sec (Mean)"] ./ sixty_initial
durationdf[!, "1Hz 60sec (SD)"] .= durationdf[!, "1Hz 60sec (SD)"] ./ sixty_initial
durationdf[!, "1Hz 90sec (Mean)"] .= durationdf[!, "1Hz 90sec (Mean)"] ./ ninety_initial
durationdf[!, "1Hz 90sec (SD)"] .= durationdf[!, "1Hz 90sec (SD)"] ./ ninety_initial

# Save the normalized data to a new CSV file
CSV.write(joinpath(datadir, "CaMKAR-duration-normalized.csv"), durationdf)

# ## Chemical data
chemicaldf = CSV.read(joinpath(datadir, "CaMKAR-chemical.csv"), DataFrame)
ctl_initial = chemicaldf[!, "Ctrl Mean"][1]
iso_initial = chemicaldf[!, "isoproterenol 100nM Mean"][1]
caf_initial = chemicaldf[!, "caffeine 20mM Mean"][1]
h2o2_50_initial = chemicaldf[!, "H2O2 50uM Mean"][1]
h2o2_100_initial = chemicaldf[!, "H2O2 100uM Mean"][1]
h2o2_200_initial = chemicaldf[!, "H2O2 200uM Mean"][1]
as10093_initial = chemicaldf[!, "AS 10093 Mean"][1]
cala_initial = chemicaldf[!, "CalA Mean"][1]

# Normalize the mean and SD relative to the initial value
chemicaldf[!, "Ctrl Mean"] .= chemicaldf[!, "Ctrl Mean"] ./ ctl_initial
chemicaldf[!, "Ctrl SD"] .= chemicaldf[!, "Ctrl SD"] ./ ctl_initial
chemicaldf[!, "isoproterenol 100nM Mean"] .= chemicaldf[!, "isoproterenol 100nM Mean"] ./ iso_initial
chemicaldf[!, "isoproterenol 100nM SD"] .= chemicaldf[!, "isoproterenol 100nM SD"] ./ iso_initial
chemicaldf[!, "caffeine 20mM Mean"] .= chemicaldf[!, "caffeine 20mM Mean"] ./ caf_initial
chemicaldf[!, "caffeine 20mM SD"] .= chemicaldf[!, "caffeine 20mM SD"] ./ caf_initial
chemicaldf[!, "H2O2 50uM Mean"] .= chemicaldf[!, "H2O2 50uM Mean"] ./ h2o2_50_initial
chemicaldf[!, "H2O2 50uM SD"] .= chemicaldf[!, "H2O2 50uM SD"] ./ h2o2_50_initial
chemicaldf[!, "H2O2 100uM Mean"] .= chemicaldf[!, "H2O2 100uM Mean"] ./ h2o2_100_initial
chemicaldf[!, "H2O2 100uM SD"] .= chemicaldf[!, "H2O2 100uM SD"] ./ h2o2_100_initial
chemicaldf[!, "H2O2 200uM Mean"] .= chemicaldf[!, "H2O2 200uM Mean"] ./ h2o2_200_initial
chemicaldf[!, "H2O2 200uM SD"] .= chemicaldf[!, "H2O2 200uM SD"] ./ h2o2_200_initial
chemicaldf[!, "AS 10093 Mean"] .= chemicaldf[!, "AS 10093 Mean"] ./ as10093_initial
chemicaldf[!, "AS 10093 SD"] .= chemicaldf[!, "AS 10093 SD"] ./ as10093_initial
chemicaldf[!, "CalA Mean"] .= chemicaldf[!, "CalA Mean"] ./ cala_initial
chemicaldf[!, "CalA SD"] .= chemicaldf[!, "CalA SD"] ./ cala_initial

# Save the normalized data to a new CSV file
CSV.write(joinpath(datadir, "CaMKAR-chemical-normalized.csv"), chemicaldf)
