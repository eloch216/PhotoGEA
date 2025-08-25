## new release

### Releasing a new version of PhotoGEA

When releasing a new version of PhotoGEA, it is essential to make sure that
`NEWS.md` is up to date and that the new version is acceptable to CRAN. This can
be achieved through the following steps:

1. Choose the next version number; our conventions for semantic versioning are
   described in `NEWS.md`.

2. Make a new release branch whose name is formatted as `release-vX.Y.X`, where
   `X.Y.Z` is the new version number, as described in
   `developer_documentation/git_branching_model.md`.

3. Update the `DESCRIPTION` file with the new release date (set to today's
   date), the new version number, and any necessary changes to the author list.

4. Change the `UNRELEASED` header in `NEWS.md` to the new version number, and
   check the contents of this section to ensure it is clear and complete. It may
   be helpful to look through the list of completed pull requests on GitHub to
   check for any important changes that may have been missed.

5. Run `R CMD check` to see any NOTEs that are reported. Check the contents of
   `cran-comments.md` to make sure it accurately reflects the `R CMD check`
   notes.

6. Make a pull request, requesting to merge the new branch into `develop`. There
   is one extra requirement for a release branch: the branch should not be
   merged until the new version has been accepted by CRAN. Sometimes CRAN may
   have issues with a new release, and it is better to address them _before_
   finalizing the release. (Otherwise, there may be several releases in quick
   succession with only minor or trivial changes between them; for example,
   BioCro versions 3.1.1 and 3.1.2.) See the "Submitting to CRAN" section below
   for instructions.

7. When the new version is on CRAN, merge the release branch into `develop` and
   then `main`.

### Submitting a new version of PhotoGEA to CRAN

The package maintainer is responsible for submitting to CRAN, and the process
consists of the following steps:

1. Make sure `R CMD check` does not produce warnings on any operating system or
   version of R, especially the current release version and the current
   development version on Linux and Windows. Our GitHub actions automatically
   run this for all pull requests.

   To be extra safe, it is also recommended to use the
   [win-builder service](https://win-builder.r-project.org/) to check the
   package again on the development version of R for Windows. The easiest way to
   do this is to call `devtools::check_win_devel()` from an R session running in
   the PhotoGEA root directory. Note that all email correspondance will be sent
   to the package maintainer regardless of who performs this action.

2. Submit the package, cross fingers, and wait. There are two ways to submit:

   1. Call `devtools::submit_cran()` from an R session running in the PhotoGEA
      root directory. This method will build the package and submit it to CRAN,
      automatically including the contents of `cran-comments.md` along with the
      built package.

   2. Submit manually. First, build the package with
      `R CMD build`. Second, attach the resulting `.tar.gz` file to the
      [CRAN form](https://cran.r-project.org/submit.html). Finally, paste the
      contents of `cran-comments.md` in the form's comments box and submit.

   With either method, the package maintainer will receive an email message
   with a link that must be followed to confirm the submission.

3. If any issues are found by CRAN, address them and try again.

   - If the checks fail only from permissible NOTES, such as using C++11, reply
     to the email indicating the justification, for example, "We use a library
     that uses C++11." You can restate what is in `cran-comments.md`.

   - If the package fails a manual review and must be resubmitted, make sure to
     increment the PhotoGEA package version, or there will be a warning about an
     insufficient package version.

PhotoGEA's online testing system should catch most issues before reaching this
point, but sometimes CRAN starts enforcing rules that are not clearly explained
anywhere or not checked by `R CMD check`. The most up-to-date description of
CRAN's requirements can be obtained from the following official and
semi-official sources, which each offer a different perspective on CRAN
submission:

 - https://cran.r-project.org/web/packages/policies.html

 - https://r-pkgs.org/release.html#release

 - https://contributor.r-project.org/cran-cookbook/code_issues.html

 - https://github.com/DavisVaughan/extrachecks
