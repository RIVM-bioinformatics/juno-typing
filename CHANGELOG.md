# Changelog

## [0.7.1](https://github.com/RIVM-bioinformatics/juno-typing/compare/v0.9.1...v0.7.1) (2025-08-04)


### Features

* add context to salmonella serotyping report ([2123f2e](https://github.com/RIVM-bioinformatics/juno-typing/commit/2123f2efa79e4bb12a3075d62fe71af9d08bcb22))
* add csv with extra hits ([ba95c91](https://github.com/RIVM-bioinformatics/juno-typing/commit/ba95c913cade2a7f3e191185e76896699395a6f5))
* add pertussis vaccine mlst ([c900f3a](https://github.com/RIVM-bioinformatics/juno-typing/commit/c900f3a3d58a43137fcab237884cdc247e369b3d))
* add sonneityping ([82afa52](https://github.com/RIVM-bioinformatics/juno-typing/commit/82afa521ce553d1dd8c451fbbe4fc1d335e0be7a))


### Bug Fixes

* Add checkout step so release-please can access latest git tag ([e8e7969](https://github.com/RIVM-bioinformatics/juno-typing/commit/e8e7969eb420609011b39fe4bfca4cf0b97e34df))
* add non-typed samples in multireport ([a212dcb](https://github.com/RIVM-bioinformatics/juno-typing/commit/a212dcb7377669d8866ad10e260f7c153f209c29))
* Added missing placeholder file ([35f62bc](https://github.com/RIVM-bioinformatics/juno-typing/commit/35f62bcb8d34efb57772410b7f59affc3beddc40))
* bump shigatyper to 2.0.5 ([92f8034](https://github.com/RIVM-bioinformatics/juno-typing/commit/92f803403915610d620024a213060b3859ed9ed9))
* catch attribute error if metadata is mssing ([59c6839](https://github.com/RIVM-bioinformatics/juno-typing/commit/59c6839bf339713f035915b23d05347ed012cd80))
* channel order ([abef539](https://github.com/RIVM-bioinformatics/juno-typing/commit/abef53977abef671f1acb2640404a78dbd3e0878))
* conda env file name ([a2aa62d](https://github.com/RIVM-bioinformatics/juno-typing/commit/a2aa62d56b9ac39feb927aef0ece3d28a7a536dc))
* Convert O_gene into a str ([02209e7](https://github.com/RIVM-bioinformatics/juno-typing/commit/02209e74cac3c9585648f286209eeea341c70135))
* Copy barrnap input temporarily ([711d186](https://github.com/RIVM-bioinformatics/juno-typing/commit/711d1862a7423bba0bd0602d911a3f8259919629))
* create expected output when shiga check fails ([afcfc16](https://github.com/RIVM-bioinformatics/juno-typing/commit/afcfc161f9342a32ac23a6d2a8083b2f4fb7ba8d))
* disable --update on grid ([9e46912](https://github.com/RIVM-bioinformatics/juno-typing/commit/9e469124eaebefbf3ce36e1e0660fc1eb3c0eb63))
* dont check for blast file ([00dd144](https://github.com/RIVM-bioinformatics/juno-typing/commit/00dd144e444fee61eefdb978ad60a8b9a5aa9530))
* fix bash if statement ([c87db73](https://github.com/RIVM-bioinformatics/juno-typing/commit/c87db73c208510940543f2e83bb5e684d2d36c7d))
* Fix channel order in master env ([3d9ec79](https://github.com/RIVM-bioinformatics/juno-typing/commit/3d9ec79679b69bc6962895afa67ea778c9adb534))
* fix check for existing mlst db ([fc250c2](https://github.com/RIVM-bioinformatics/juno-typing/commit/fc250c2f7e8b880a1bfeb14f45d302afcfd83af3))
* fix DAG logic for non shigella typings ([3aa268c](https://github.com/RIVM-bioinformatics/juno-typing/commit/3aa268ca8b4a2c6670e7527ed8ff051b28b15317))
* fix error when no shigella samples observed ([b2a5538](https://github.com/RIVM-bioinformatics/juno-typing/commit/b2a55389109613a6f98145958f84edddfbbae00f))
* Fix libarchive version ([463557c](https://github.com/RIVM-bioinformatics/juno-typing/commit/463557c4901f360603cd63f01021cdb982a2814f))
* fix pulp version ([6fb2eac](https://github.com/RIVM-bioinformatics/juno-typing/commit/6fb2eaccd2d086c2490eb69655f679bd1af21477))
* fix seroba db rule ([3c47dd4](https://github.com/RIVM-bioinformatics/juno-typing/commit/3c47dd444032564d0a405da7a6fd5c3f8f3e8b66))
* Fix Seroboba untyped output breaking multireport ([29edbf7](https://github.com/RIVM-bioinformatics/juno-typing/commit/29edbf756a6ed8247381e5af47fc006469a19f69))
* grep correct data from stdout ([7811a4f](https://github.com/RIVM-bioinformatics/juno-typing/commit/7811a4f5bf2c06434203144cc14f312772d134f6))
* path to bordetella db ([add37d3](https://github.com/RIVM-bioinformatics/juno-typing/commit/add37d3c485e5e20585dbd9b661d30cbd80e14c4))
* remove biocore channel and replace blast-plus with blast ([40625dd](https://github.com/RIVM-bioinformatics/juno-typing/commit/40625dd1c904169da2c9c9eff1a350e7a43065c9))
* remove intel channel ([955a864](https://github.com/RIVM-bioinformatics/juno-typing/commit/955a864a032c4a9e6a6b00d180dd94a12c07e8ba))
* replaced cconcisus_curvus with cconcisus ([3732ea2](https://github.com/RIVM-bioinformatics/juno-typing/commit/3732ea28fb79b1fbc71daa093af5e5dd341e00f6))
* revert old check ([c2d6019](https://github.com/RIVM-bioinformatics/juno-typing/commit/c2d601928bc22f634c0f2ca90ed920e315664da2))
* run checkpoint locally ([15edca0](https://github.com/RIVM-bioinformatics/juno-typing/commit/15edca0f6fe10a3e0c0cdeae85e8765d71f90c8a))
* set bordetella scheme name in object ([4b78d67](https://github.com/RIVM-bioinformatics/juno-typing/commit/4b78d6746ee24058350f706864858ab6cb5bd491))
* Set channel priority to strit in run_pipeline ([a72253b](https://github.com/RIVM-bioinformatics/juno-typing/commit/a72253be6735f9afacdcdddc1bb683d7fc699c5b))
* set correct db check seroba ([3a83e08](https://github.com/RIVM-bioinformatics/juno-typing/commit/3a83e08aad4a2f197811cb9548c2b83f546db7c2))
* set sample-specific skeleton dir ([be75ef7](https://github.com/RIVM-bioinformatics/juno-typing/commit/be75ef7c35258d72e2fc9d36d0baaf1788afb28b))
* stop downloading neisseria database, neisseria database is now locally available ([2c0fcea](https://github.com/RIVM-bioinformatics/juno-typing/commit/2c0fcea6416129391a546ff2fac4488335c0c780))
* update c. acnes mlst scheme name ([bee5f08](https://github.com/RIVM-bioinformatics/juno-typing/commit/bee5f08105dadcdd84d494ac8ffda9fefde0d66e))
* Update juno library to v2.2.2 ([bae06c0](https://github.com/RIVM-bioinformatics/juno-typing/commit/bae06c018448f644154990fed1f85f0d8ddf18ad))
* Update master_env to fix mamba crash ([293b9ed](https://github.com/RIVM-bioinformatics/juno-typing/commit/293b9ed2bf16f5493a150bc31d67e4b0548cb238))
* Update SeqSero context notes ([244936f](https://github.com/RIVM-bioinformatics/juno-typing/commit/244936fc5002b30720b9e0bee390fef083645232))
* update shigatyper ([828fbe5](https://github.com/RIVM-bioinformatics/juno-typing/commit/828fbe5f70d03d83f9d488229253e6fdcc04ba0e))
* use same mamba version in CI as pipeline ([72d88b6](https://github.com/RIVM-bioinformatics/juno-typing/commit/72d88b6ab2af6e16aef32f9597ab64ab3d17baed))
* working version locally ([e8c3b59](https://github.com/RIVM-bioinformatics/juno-typing/commit/e8c3b592b38655555c866c6a7d5449c3e3a207de))


### Dependencies

* Add explicit libarchive version ([f4ab775](https://github.com/RIVM-bioinformatics/juno-typing/commit/f4ab7755cc8f73b9943781fb9feeb1cb29b6c449))
* fix seqsero install ([baf735b](https://github.com/RIVM-bioinformatics/juno-typing/commit/baf735b4c1f6b43b4a4acf2077302a319fadc578))
* remove anaconda and defaults channel and add nodefaults channel, also update dependencies in multiple envs ([ac6831c](https://github.com/RIVM-bioinformatics/juno-typing/commit/ac6831cc153f8333aa14f868b382d850c312d74e))
* Remove anaconda channel from environments ([3447b49](https://github.com/RIVM-bioinformatics/juno-typing/commit/3447b49faa46f350d9256bf4eaa561edba278a7e))
* Update dependencies ([9fe2b94](https://github.com/RIVM-bioinformatics/juno-typing/commit/9fe2b949018d18a9b43e349d1b85f32f8719ba71))


### Reverts

* Remove reverted git tag and changes (from 0.9.1 to 0.8.6) ([d32c2dd](https://github.com/RIVM-bioinformatics/juno-typing/commit/d32c2ddce55232cc94adfc7844ea962e058a7f86))


### Documentation

* Add 16S extraction documentation ([c13e4b1](https://github.com/RIVM-bioinformatics/juno-typing/commit/c13e4b168910b27911edc29426e94ac2af465be7))


### Miscellaneous Chores

* release 0.7.0 ([23155c8](https://github.com/RIVM-bioinformatics/juno-typing/commit/23155c84a4dd2be91bc6c872c615c41a7b0be102))
* release 0.7.1 ([88dbdcd](https://github.com/RIVM-bioinformatics/juno-typing/commit/88dbdcd0967ba1554bc89014e1d4890a3e7d691b))

## [0.8.4](https://github.com/RIVM-bioinformatics/juno-typing/compare/v0.8.3...v0.8.4) (2024-09-04)


### Bug Fixes

* update c. acnes mlst scheme name ([bee5f08](https://github.com/RIVM-bioinformatics/juno-typing/commit/bee5f08105dadcdd84d494ac8ffda9fefde0d66e))

## [0.8.3](https://github.com/RIVM-bioinformatics/juno-typing/compare/v0.8.2...v0.8.3) (2024-07-03)


### Bug Fixes

* remove intel channel ([955a864](https://github.com/RIVM-bioinformatics/juno-typing/commit/955a864a032c4a9e6a6b00d180dd94a12c07e8ba))

## [0.8.2](https://github.com/RIVM-bioinformatics/juno-typing/compare/v0.8.1...v0.8.2) (2024-05-21)


### Bug Fixes

* fix pulp version ([6fb2eac](https://github.com/RIVM-bioinformatics/juno-typing/commit/6fb2eaccd2d086c2490eb69655f679bd1af21477))

## [0.8.1](https://github.com/RIVM-bioinformatics/juno-typing/compare/v0.8.0...v0.8.1) (2024-04-16)


### Bug Fixes

* fix seroba db rule ([3c47dd4](https://github.com/RIVM-bioinformatics/juno-typing/commit/3c47dd444032564d0a405da7a6fd5c3f8f3e8b66))
* set correct db check seroba ([3a83e08](https://github.com/RIVM-bioinformatics/juno-typing/commit/3a83e08aad4a2f197811cb9548c2b83f546db7c2))

## [0.8.0](https://github.com/RIVM-bioinformatics/juno-typing/compare/v0.7.1...v0.8.0) (2024-04-12)


### Features

* add csv with extra hits ([ba95c91](https://github.com/RIVM-bioinformatics/juno-typing/commit/ba95c913cade2a7f3e191185e76896699395a6f5))

## [0.7.1](https://github.com/RIVM-bioinformatics/juno-typing/compare/v0.7.0...v0.7.1) (2024-02-29)


### Miscellaneous Chores

* release 0.7.1 ([88dbdcd](https://github.com/RIVM-bioinformatics/juno-typing/commit/88dbdcd0967ba1554bc89014e1d4890a3e7d691b))

## [0.7.0](https://github.com/RIVM-bioinformatics/juno-typing/compare/v0.6.3...v0.7.0) (2024-02-28)


### Miscellaneous Chores

* release 0.7.0 ([23155c8](https://github.com/RIVM-bioinformatics/juno-typing/commit/23155c84a4dd2be91bc6c872c615c41a7b0be102))

## [0.6.3](https://github.com/RIVM-bioinformatics/juno-typing/compare/v0.6.2...v0.6.3) (2024-02-22)


### Dependencies

* fix seqsero install ([baf735b](https://github.com/RIVM-bioinformatics/juno-typing/commit/baf735b4c1f6b43b4a4acf2077302a319fadc578))

## [0.6.2](https://github.com/RIVM-bioinformatics/juno-typing/compare/v0.6.1...v0.6.2) (2024-02-15)


### Bug Fixes

* disable --update on grid ([9e46912](https://github.com/RIVM-bioinformatics/juno-typing/commit/9e469124eaebefbf3ce36e1e0660fc1eb3c0eb63))
* fix check for existing mlst db ([fc250c2](https://github.com/RIVM-bioinformatics/juno-typing/commit/fc250c2f7e8b880a1bfeb14f45d302afcfd83af3))
* revert old check ([c2d6019](https://github.com/RIVM-bioinformatics/juno-typing/commit/c2d601928bc22f634c0f2ca90ed920e315664da2))

## [0.6.1](https://github.com/RIVM-bioinformatics/juno-typing/compare/v0.6.0...v0.6.1) (2024-01-17)


### Bug Fixes

* dont check for blast file ([00dd144](https://github.com/RIVM-bioinformatics/juno-typing/commit/00dd144e444fee61eefdb978ad60a8b9a5aa9530))

## [0.6.0](https://github.com/RIVM-bioinformatics/juno-typing/compare/v0.5.12...v0.6.0) (2023-11-14)


### Features

* add pertussis vaccine mlst ([c900f3a](https://github.com/RIVM-bioinformatics/juno-typing/commit/c900f3a3d58a43137fcab237884cdc247e369b3d))


### Bug Fixes

* catch attribute error if metadata is mssing ([59c6839](https://github.com/RIVM-bioinformatics/juno-typing/commit/59c6839bf339713f035915b23d05347ed012cd80))
* conda env file name ([a2aa62d](https://github.com/RIVM-bioinformatics/juno-typing/commit/a2aa62d56b9ac39feb927aef0ece3d28a7a536dc))
* path to bordetella db ([add37d3](https://github.com/RIVM-bioinformatics/juno-typing/commit/add37d3c485e5e20585dbd9b661d30cbd80e14c4))
* set bordetella scheme name in object ([4b78d67](https://github.com/RIVM-bioinformatics/juno-typing/commit/4b78d6746ee24058350f706864858ab6cb5bd491))

## [0.5.12](https://github.com/RIVM-bioinformatics/juno-typing/compare/v0.5.11...v0.5.12) (2023-07-18)


### Dependencies

* remove anaconda and defaults channel and add nodefaults channel, also update dependencies in multiple envs ([ac6831c](https://github.com/RIVM-bioinformatics/juno-typing/commit/ac6831cc153f8333aa14f868b382d850c312d74e))
* Remove anaconda channel from environments ([3447b49](https://github.com/RIVM-bioinformatics/juno-typing/commit/3447b49faa46f350d9256bf4eaa561edba278a7e))

## [0.5.11](https://github.com/RIVM-bioinformatics/juno-typing/compare/v0.5.10...v0.5.11) (2023-06-28)


### Bug Fixes

* update shigatyper ([828fbe5](https://github.com/RIVM-bioinformatics/juno-typing/commit/828fbe5f70d03d83f9d488229253e6fdcc04ba0e))

## [0.5.8](https://github.com/RIVM-bioinformatics/juno-typing/compare/v0.5.7...v0.5.8) (2023-05-26)


### Bug Fixes

* channel order ([abef539](https://github.com/RIVM-bioinformatics/juno-typing/commit/abef53977abef671f1acb2640404a78dbd3e0878))

## [0.5.7](https://github.com/RIVM-bioinformatics/juno-typing/compare/v0.5.6...v0.5.7) (2023-05-23)


### Bug Fixes

* Set channel priority to strit in run_pipeline ([a72253b](https://github.com/RIVM-bioinformatics/juno-typing/commit/a72253be6735f9afacdcdddc1bb683d7fc699c5b))
* stop downloading neisseria database, neisseria database is now locally available ([2c0fcea](https://github.com/RIVM-bioinformatics/juno-typing/commit/2c0fcea6416129391a546ff2fac4488335c0c780))

## [0.5.6](https://github.com/RIVM-bioinformatics/juno-typing/compare/v0.5.5...v0.5.6) (2023-04-26)


### Dependencies

* Add explicit libarchive version ([f4ab775](https://github.com/RIVM-bioinformatics/juno-typing/commit/f4ab7755cc8f73b9943781fb9feeb1cb29b6c449))

## [0.5.5](https://github.com/RIVM-bioinformatics/juno-typing/compare/v0.5.4...v0.5.5) (2023-04-25)


### Bug Fixes

* Copy barrnap input temporarily ([711d186](https://github.com/RIVM-bioinformatics/juno-typing/commit/711d1862a7423bba0bd0602d911a3f8259919629))
