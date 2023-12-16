## git branching model

### git-flow

* In general, PhotoGEA development follows the
  [git-flow](https://nvie.com/posts/a-successful-git-branching-model/) branching
  model, where there are two permanent branches (`main` and `unreleased`) and
  three types of temporary branches (hotfixes, features, and releases).

* For most contributors, it is only necessary to know that most changes should
  be accomplished through feature branches, which are branched off `unreleased`
  and merged back into `unreleased` via a pull request that must be approved by
  one or more other developers. The remainder of this section discusses our
  system in more detail.

* Beyond the basic description of git-flow, we have a few additional rules
  and clarifications specific to PhotoGEA development:

  - Any merge into `main` or `unreleased` must be done via a pull request. On
    the `main` branch, this requirement is enforced with GitHub branch
    protection rules and cannot be bypassed.

  - All pull requests into `main` and `unreleased` require approval before
    merging. This requirement is not enforced using GitHub branch protection
    rules, but please do not unilaterally merge a pull request without any form
    of approval from another developer. (There are two exceptions to this rule
    related to hotfix and release branches, described below.)

  - The PhotoGEA R package and repository follow semantic versioning of the form
    `major.minor.patch`. Hotfix branches should typically only increment the
    third (`patch`) number. Most release branches will increment the second
    (`minor`) number, and very rarely the first (`major`) number. Feature
    branches do not directly change the version number.

  - Hotfix and release branch names should be formatted like `type-version`,
    where `type` is either `hotfix` or `release`, and version is the new version
    number. For example, `hotfix-v3.0.3` would be a hotfix branch that
    changes the version number to `3.0.3`. Feature branch names should reflect
    their purpose, and can be anything other than `main`, `unreleased`,
    `hotfix-*`, or `release-*`.

  - Whenever the package version changes, a description of the related changes
    should be added to the changelog in `NEWS.md` as a new section titled with
    the new version number. Such version changes will happen via release and
    hotfix branches. Release branches will generally include changes made during
    the course of several feature branches, so they may require a significant
    amount of documentation in `NEWS.md`. To make it easier to prepare a list of
    these changes, each feature branch should include a description of its own
    changes in an `UNRELEASED` section of `NEWS.md`. For more information about
    updating the changelog, please see the comment at the top of `NEWS.md`.

* The following is a short description of PhotoGEA's implementation of the
  git-flow branching model:

  - The `main` branch always contains the latest stable release of the PhotoGEA
    GitHub repository. In fact, a new tag and release is made any time changes
    are merged into `main`, and such changes should always be accompanied by an
    increment to the PhotoGEA R package version number.

  - The `unreleased` branch contains bug fixes and new features that have been
    completed but may not have been released yet.

  - Hotfix branches are for urgent changes that warrant immediate incorporation
    into `main`; typically these are bug fixes. A hotfix branch should be
    branched off `main`. When ready, it should first merged into `main` via an
    approved pull request. Then, a second pull request should be made to merge
    into `unreleased`. If there are no merge conflicts or test failures, this
    second request can be merged without any additional approval. When both
    merges are complete, the hotfix branch should be deleted.

  - Feature branches are for less urgent changes that do not require immediate
    release; for example, the addition of a new R function to the package
    namespace. A feature branch should be branched off `unreleased` and merged
    back into `unreleased` via an approved pull request. Before merging,
    `NEWS.md` should be updated with a description of the changes under the
    `# UNRELEASED` section (see `NEWS.md` for more details). When the merge into
    `unreleased` is complete, the feature branch should be deleted.

  - When sufficiently many changes have accumulated in `unreleased` to justify a
    new version of the package, a release branch is used to move the changes
    from `unreleased` to `main`. Release branches should not include substantive
    changes; rather, a release branch is primarily used to increment the package
    version and to ensure an up-to-date changelog in `NEWS.md`. A release branch
    should be branched off `unreleased`. When it's ready, it should first merged
    into `main` via an approved pull request. Then, a second pull request should
    be made to merge into `unreleased`. If there are no merge conflicts or test
    failures, this second request can be merged without any additional approval.
    When both merges are complete, the release branch should be deleted.
