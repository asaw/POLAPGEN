#!/bin/bash
qsub mqsub_1_ref_cow_grid.sh
qsub -W depend=afterok:xxxxxxxxx.ce.reef mqsub_2_ref_cow_opt.sh
qsub -W depend=afterok:xxxxxxxxx.ce.reef mqsub_3_cow_grid.sh
qsub -W depend=afterok:xxxxxxxxx.ce.reef mqsub_4_cow_opt.sh
qsub -W depend=afterok:xxxxxxxxx.ce.reef run_5_table.sh
qsub -W depend=afterok:xxxxxxxxx.ce.reef run_6_table_peak_division.sh