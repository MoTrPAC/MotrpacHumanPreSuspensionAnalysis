# Retired: precovid-repro builds this now.
#
#   precovid-repro/scripts/10_build_data/13_build_fcm/FCM_clustering_results.R
#
# The objects this script used to write into data/ are produced by the
# pipeline and carried in by its Stage 3 (30_update_relevant_packages), so
# running anything here would write a second, unversioned copy from inputs
# the release no longer reads. The code is removed rather than left runnable
# beside the thing that replaced it: two live generators for one object is
# how a package ends up shipping data nobody can trace.
#
# inst/PROVENANCE.tsv records, per object, which pipeline step built it and
# whether it was regenerated or staged verbatim.

