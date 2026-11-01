# MicroCyte v1.0.0

MicroCyte is a guided microcytometry workflow developed by Jason Needham at the University of
Alabama at Birmingham (2020–2023). It organizes microscopy images, prepares them for analysis,
runs cell and region-of-interest measurements, and combines the results into analysis-ready data.

![MicroCyte workflow](https://github.com/SThompsonLab/MicroCyte/blob/main/MicroCyte.png)

For a complete worked example, see
`MicroCyte_Example_Walkthrough_for_Cell_Cycle.pdf` in the MicroCyte folder.

## What you need

MicroCyte currently requires three components:

1. **R** — download it from [The Comprehensive R Archive Network](https://cran.r-project.org/).
2. **RStudio Desktop** — optional, but recommended for an easier experience. Download it from
   [Posit](https://posit.co/download/rstudio-desktop/).
3. **Fiji/ImageJ** — required for detecting cells and measuring the microscopy images. Download it
   from [fiji.sc](https://fiji.sc/).

MicroCyte works with the same R files on Windows and macOS. Keep the entire downloaded MicroCyte
folder together; do not move `microCyte.R`, `schema.csv`, or the `bin` folder away from one another.

## Download MicroCyte

1. On the MicroCyte GitHub page, select **Code**, then **Download ZIP**.
2. Extract the ZIP file to a location where you have permission to create and modify files.
3. Open the extracted `MicroCyte` folder.

Do not run MicroCyte from inside the ZIP file.

## Start MicroCyte

1. Open R or RStudio.
2. Set the working directory to the extracted MicroCyte folder—the folder containing
   `microCyte.R` and `schema.csv`.
   - In RStudio, open `microCyte.R`, then select
     **Session → Set Working Directory → To Source File Location**.
3. Install the MicroCyte analysis package from the project directory:

```r
install.packages(".", repos = NULL, type = "source")
```

4. Start the workflow CLI:

```r
source("microCyte.R")
run_microcyte()
```

The MicroCyte menu will appear. Installing the package makes the analysis functions available to
menu option **9 — analyze** and through `library(MicroCyte)`.

### Install the R dependencies

The first time you run MicroCyte, choose menu option **0 — Setup**. MicroCyte will check for its
required R packages and download any that are missing. An internet connection is required during
this step. Package installation may take several minutes.

After setup finishes, return to the menu and continue with the workflow. Setup normally only needs
to be run again after reinstalling or updating R.

## Typical workflow

| Order | Menu option | Purpose |
|---:|---|---|
| 1 | **0 — Setup** | Install missing R packages. |
| 2 | Edit `schema.csv` | Describe the samples, image channels, and experimental metadata. |
| 3 | **1 — dir_gen** | Create the experiment and image directories. |
| 4 | Add images | Place the microscopy images in the generated directories. |
| 5 | **2 — purgo** | Prepare images and optionally correct spectral overlap. |
| 6 | **3 — re_name** | Rename image channels using `schema.csv`. |
| 7 | Fiji/ImageJ macro | Detect cells and export measurements. |
| 8 | **4 — ima_gen** | Combine the ImageJ measurements for each image. |
| 9 | **5 — concato** | Combine images belonging to the same condition. |
| 10 | **6 — unite** | Combine conditions into final experiment datasets. |
| 11 | **9 — analyze** | Load the included analysis functions. |

Menu option **7 — full run** guides you through the complete sequence. The workflow pauses when you
need to collect images or run the Fiji/ImageJ macro.

## Prepare `schema.csv`

The schema describes the experiment. Each row represents one unique sample or condition.

### Required information

- `notebook_id` must uniquely identify each row.
- Channel columns such as `CH_1`, `CH_2`, and `CH_3` map the microscope's channel labels to the
  biological target names.
- A DNA or nuclear stain such as DAPI must use the target name `dna`. The ImageJ macro uses this
  channel to identify nuclei.
- Additional columns may contain metadata such as time point, infection status, treatment, or
  biological replicate.

Example:

| notebook_id | CH_1 | CH_2 | CH_3 | CH_4 | time | infection |
|---|---|---|---|---|---|---|
| JNp1_1 | betaActin | lamp1 | dna | nucleocapsid | 24hpi | infected |
| JNp1_2 | betaActin | lamp1 | dna | nucleocapsid | 48hpi | infected |
| JNp1_3 | betaActin | lamp1 | dna | nucleocapsid | 48hpi | uninfected |

When generating directories, MicroCyte can use either:

- **number** — creates directories from `notebook_id`, such as `JNp1_1`; or
- **name** — creates descriptive directories from the metadata, such as
  `48hpi_uninfected`.

Every row must have a unique `notebook_id` or a unique combination of metadata values.

## Generate the image directories

Choose menu option **1 — dir_gen**.

MicroCyte creates the following directories if they do not already exist:

- `files` — source images, intermediate files, and condition-level results;
- `data` — final combined datasets; and
- `figures` — generated plots and figures.

You will be asked whether images are captured automatically or individually. Choose **single** to
have MicroCyte create numbered image folders for you. Choose **auto** if the microscope already
creates a separate folder for each field or image set.

## Add the microscopy images

Place each image set in its generated directory under `files`.

- Use TIF images.
- Save each fluorescent channel as a separate image rather than a merged image.
- Images may use the microscope channel names, such as `1_CH_1.tif`, or their target names.
- Keep all channels from one field of view in the same image directory.
- Ensure the nuclear channel is mapped to `dna` in `schema.csv`.

## Prepare and rename the images

Choose menu option **2 — purgo** to convert the images and configure spectral-overlap correction.
Available modes are:

- **none** — prepare the file structure without spectral correction;
- **auto** — use the default automatic correction;
- **manual** — choose the overlapping channels; and
- **full** — apply sequential correction across the configured channels.

MicroCyte can generate quality-control plots showing relationships between channel intensities.
After Purgo finishes, it offers to run the rename step automatically.

Choose menu option **3 — re_name** separately if you need to rename the images again after editing
`schema.csv`.

## Run Fiji/ImageJ

The Fiji macro is located at:

`bin/scripts/scripts/4_YggData_multi.ijm`

Open the macro in Fiji and run it after the images have been prepared and renamed. The macro:

1. Uses the `dna` image to identify nuclei.
2. Measures every fluorescent channel within the nuclear regions.
3. Optionally measures non-nuclear and perinuclear regions.
4. Saves CSV measurement files into each image directory.

Test the detection settings on several representative images before processing the full experiment.
The nuclear intensity threshold, size range, circularity, and target-specific settings are defined at
the beginning of the macro.

## Combine the results

After the Fiji macro finishes, return to MicroCyte:

1. Choose **4 — ima_gen** to associate measured regions with nuclei and create per-image datasets.
2. Choose **5 — concato** to combine the images within each condition.
3. Choose **6 — unite** to combine conditions into final experiment datasets.

During Concato, MicroCyte can normalize measurements between image sets. During Unite, you can use
all detected cells or draw the same number of cells from each condition. Conditions can be combined
all together or grouped by a value from the schema.

Final combined CSV files are written to the `data` directory.

## Analysis tools

The analysis tools are provided by the `MicroCyte` R package. Load them directly with:

```r
library(MicroCyte)
```

Menu option **9 — analyze** loads the same exported functions into the interactive analysis
environment for compatibility with the original workflow.

Analysis functions now require their data and paths explicitly and return their results rather than
creating or replacing objects in the global environment. For example:

```r
cells <- read.csv("data/experiment.csv")
plaque_result <- foci_caller(cells)
summary <- explore(cells, scheme = "schema.csv", save_file = FALSE)
density_result <- color_density(cells, x = "dna_norm", y = "edu_norm")
```

The package includes:

- **SirMixaPlot tools** — explicit data transformations and returned ggplot objects;
- **CAP-IFA tools** — identify and summarize infected plaques;
- **explore** and **plumeria** — create schema-based summaries and normalized metrics;
- **icellate**, heat-map, and correlation tools — process images using explicit paths;
- **ternaplot**, **ergodic_slicing**, and **smoke_stack** — return plots and transformed data; and
- **pull_n_zip** — collect and archive raw images without changing the working directory.

Specialized dependencies are optional and declared in `DESCRIPTION`. A function reports its missing
package when that capability is called, so an unavailable image or ternary plotting dependency no
longer prevents unrelated analysis tools from loading.

The final CSV files can also be opened in spreadsheet software, R, or flow-cytometry software that
supports CSV input.

## If something goes wrong

- Confirm that R's working directory is the main MicroCyte folder.
- Confirm that `schema.csv`, `microCyte.R`, and the `bin` folder are still together.
- Confirm that every `notebook_id` is unique.
- Confirm that the nuclear target is named `dna`.
- Confirm that each image set contains the expected TIF channel images.
- Run menu option **0 — Setup** again if R reports that a package is missing.
- Check the quality-control images before trusting the final measurements.

MicroCyte is intended to make high-content microscopy analysis easier and more reproducible.
