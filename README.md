![metahit-removebg-preview](https://github.com/user-attachments/assets/205507ac-2766-470e-9c6d-2ddebc279f74)


## Software Setup
Run `bash setup.sh` to install external dependenceies.
Run `conda env create -f env.yaml` in root directory
Run `conda activate metahit_env` to activate conda env.
RUN `conda env remove -n metahit_env` to delete env
Then, in the conda environment, directly run `bash demo.sh` and the pipeline will be triggered.


## Database setup

The folder (`db_setup`) contains three scripts to download and set up databases for **CheckM**, **CheckM2**, and **GTDB-Tk**. By default, each script downloads the database into a `database` folder in your current working directory. You can optionally provide a custom path.

---


- **CheckM database downloading:** `checkm_db.sh`
  Usage:
  ```bash
  bash db_setup/checkm_db.sh [DB_DIR]
  ```

- **CheckM2 database downloading:** `checkm2_db.sh`
  Usage
  ```bash
  bash db_setup/checkm2_db.sh [DB_DIR]
  ```

- **GTDB-TK database downloading:** `tdbtk_db.sh`
  Usage:
  ```bash
  bash db_setup/gtdbtk_db.sh [DB_DIR] [RELEASE_VERSION]
  ```
  Example:
  ```bash
  bash db_setup/gtdbtk_db.sh /your/custom/path 214
  ```