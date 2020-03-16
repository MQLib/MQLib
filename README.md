![MQLib](mqlib.png)

# MQLib

* An implementation of dozens of heuristics for the **M**ax-cut and **Q**UBO combinatorial optimization problems.
* A machine learning-based hyper-heuristic that tries to select the best heuristic for a given instance.
* Scripts to evaluate heuristics on Amazon EC2 and to analyze the results.

This library and the related systematic heuristic evaluation strategy are described in [the paper](paper/SECM_final.pdf). To cite the MQLib, please use:
```
@article{DunningEtAl2018,
  title={What Works Best When? A Systematic Evaluation of Heuristics for Max-Cut and {QUBO}},
  author={Dunning, Iain and Gupta, Swati and Silberholz, John},
  year={2018},
  journal={{INFORMS} Journal on Computing},
  volume={30},
  number={3}
}
```

This code is licensed under the MIT License - see `/LICENSE` for details. MQLib was created by [Iain Dunning](http://iaindunning.com), [Swati Gupta](https://swatigupta.tech), and [John Silberholz](http://johnsilberholz.com).

## Obtaining Source Code and Building the MQLib

The first step to building the MQLib is to download the source code by cloning the GitHub repository. This can be done by running `git clone https://github.com/MQLib/MQLib` from the command line.

After changing directories to the MQLib directory, run the `make` command to build the project. Currently, this requires several GNU compiler tools (the `g++` compiler and the `ar` archive utility) as well as several common unix utilities (`find`, `sed`, `rm`, and `make`). These tools should be readily available on *nix distributions as well as on Mac operating systems with the XCode tools installed. However, they can also be accessed on Windows machines through Cygwin.

After this compilation step, an executable should be available at `bin/MQLib` and a linkable library should be available at `bin/MQLib.a`.

## Use Cases

The MQLib can be used for a number of different purposes. We provide pointers for a number such use cases here:

### Running MQLib Heuristics by Executable or by Linking to the Library

The MQLib provides access to a state-of-the-art hyper-heuristic as well as a number of heuristics from the literature. These heuristics can be run either from the command line or by linking to the MQLib. Details are available in the [Using the MQLib guide](bin/README.md).

### Contributing a New Heuristic

A central goal of the MQLib is to provide open-source implementations of many Max-Cut and QUBO heuristics. As a result, we are excited to add new heuristics and their results into the library. To do so, use the following steps:

1. Create a fork of the MQLib repository, as described in the "Contributing to the MQLib" section of this document.
2. Implement your new heuristic in C++. Details of how to implement a new heuristic are available in the [developer's guide](src/README.md). The heuristic can be tested using the command-line interface to the MQLib, as detailed in the [Using the MQLib guide](bin/README.md). It may be useful to test the new heuristic on problem instances with a range of properties, which can be identified in [data/metrics.csv](data/metrics.csv) and downloaded with [scripts/downloadGraph.py](scripts/downloadGraph.py).
3. Commit your new code to your fork of the MQLib repository, as described in the "Contributing to the MQLib" section of this document. If desired, the MQLib team would be glad to review any code at this stage and to provide suggestions before extensive testing on the cloud. Please email Iain Dunning (iaindunning@gmail.com), Swati Gupta (swatig@gatech.edu), and John Silberholz (john.silberholz@gmail.com) with any such requests.
4. Test your new heuristic using Amazon Web Services, as detailed in the [Reproducible Parallel Computation with Amazon Web Services guide](Cloud/README.md). Make sure to specify your forked repository when doing that run and to change the `data/heuristics.txt` file to test your new heuristic only.
5. You should uncomment the appropriate lines under "Evaluation results" in [analysis/analyze.R](analysis/analyze.R) to combine your new results with the existing results. Adding new results may affect multiple components of the analysis. The largest initial impact is that "interesting" subset will change. Naturally this will also affect Table 3 (summary of results) and Figure 7 (plots of heuristic performance on standard and interesting). Finally, the hyper-heuristic itself will change.
6. Use `git add` to add to the repository a zipped version of the heuristic results as a new file in the `data` folder and the new hyper-heuristic data file corresponding to your new heuristic in the `hhdata` folder. As described in the "Contributing to the MQLib" section of this document, commit these new files as well as an updated version of the [`data` folder readme file](data/README.md) describing the computational experiments, an updated version of [analysis/analyze.R](analysis/analyze.R) modified to load your new computational results, and the updated hyper-heuristic files in the `hhdata` folder.
7. As described in the "Contributing to the MQLib" section of this document, submit a pull request with all committed modifications to the MQLib repository.

### Contributing a New Problem Instance Metric

The MQLib uses fast-to-compute problem instance metrics to predict the best-performing heuristic or set of heuristics for a given problem instance. This prediction task forms the core of the hyper-heuristic shipped with the MQLib. We welcome any new fast-to-compute metrics that can improve the performance of the hyper-heuristic; to keep the hyper-heuristic speedy, new metrics should not be much slower than the slowest metric currently used in the project. To add a new metric or set of metrics to the project, use the following steps:

1. Create a fork of the MQLib repository, as described in the "Contributing to the MQLib" section of this document.
2. Implement your new metric (or set of metrics) in C++. Details of how to implement a new metric are available in the [developer's guide](src/README.md). The metric can be tested using the command-line interface to the MQLib, as detailed in the [Using the MQLib guide](bin/README.md). It may be useful to test the metric computation code on problem instances with a range of properties, which can be identified in [data/metrics.csv](data/metrics.csv) and downloaded with [scripts/downloadGraph.py](scripts/downloadGraph.py).
3. Commit your new code to your fork of the MQLib repository, as described in the "Contributing to the MQLib" section of this document.
4. Compute your new metric for all problem instances using Amazon Web Services, as detailed in the [Reproducible Parallel Computation with Amazon Web Services guide](Cloud/README.md). Make sure to specify your forked repository when doing that run.
5. Follow the instructions in [analysis/analyze.R](analysis/analyze.R) for combining your new metric with the exisiting metrics, if you have not already generated a new file containing both the old metrics and your new metrics for all instances.
6. Run through the analysis with your new metric. Relevant sections include the performance of the hyperheuristic with the benefit of your new metric, and the interpretable models.
7. As described in the "Contributing to the MQLib" section of this document, commit an updated version of [data/metrics.csv](data/metrics.csv) containing the new metrics, an updated version of the [`data` folder readme file](data/README.md) describing the updated metrics, and the updated hyper-heuristic files in the `hhdata` folder.
8. As described in the "Contributing to the MQLib" section of this document, submit a pull request with all committed modifications to the MQLib repository.

### Other use cases

There are other modifications possible for the MQLib that would require more involved changes to the code base or related `mqlibinstances` S3 bucket (which stores the problem instances). We summarize them here and encourage researchers to reach out to the MQLib team before working to implement these sorts of extensions by emailing Iain Dunning (iaindunning@gmail.com), Swati Gupta (swatig@gatech.edu), and John Silberholz (john.silberholz@gmail.com).

* **Extending the set of problem instances used for testing:** The MQLib team is enthusiastic about extending the set of problem instances used in testing Max-Cut and QUBO heuristics, and especially in adding problem instances from real-world applications of the problems. Because the `mqlibinstances` S3 bucket used for sharing instances is not publicly writable, researchers should contact the MQLib team to add new instances.
* **Testing across heuristic parameter settings:** Such a modification would require significant changes to the constructors of all heuristics to enable parameters to be specified. Further changes would be needed to [src/main.cpp](src/main.cpp) to enable command-line specification of parameters and to [Cloud/MQLibDispatcher.py](Cloud/MQLibDispatcher.py) and [Cloud/MQLibRunner.py](Cloud/MQLibRunner.py) to specify parameters during execution on the cloud.

## Contributing to the MQLib

As detailed above, the first step in contributing to the MQLib involves forking the GitHub repository. You can do this using the following steps:

* Log into [github.com](https://github.com) and go to the project page, [https://github.com/MQLib/MQLib](https://github.com/MQLib/MQLib).
* Click the "Fork" button, selecting yourself.
* Navigate to your forked repository by clicking your name at the top right of the screen, clicking the "Repositories" tab, and then clicking the MQLib link.
* Select the URL in the "HTTPS clone URL" section on the right side of the screen, and copy it to your clipboard.
* On the command line, navigate to where you want to store the repository and run `git clone URL`, where `URL` is the URL that you copied in the previous step.

To commit changes to your fork (which is needed if you have made any changes to the source code of the project and want those changes to be in effect when performing a parallel run with Amazon Web Services), you can use the following steps:

* Add new files with the `git add` command; for instance, if you implemented a new heuristic you would add the header and source files with commands like `git add include/heuristics/maxcut/silberholz2015.h` and `git add src/heuristics/maxcut/silberholz2015.cpp`.
* Commit your changes with the command `git commit -a -m "MESSAGE"`, leaving a descriptive message.
* Push your edits with the command `git push`.

Note that you may need to commit to your fork several times. For instance, when implementing a new heuristic you would need to commit before cloud-based testing (to add the heuristic to the repository) and then after analysis is complete (to share computational results and the updated hyper-heuristic). When implementing a new metric you would need to commit before cloud-based metric computation (to add the metric source code to the repository) and then after analysis is complete (to share the updated metric values and the updated hyper-heuristic).

When you are ready to integrate your fork into the main repository to share your code and results with others, please use the following steps:

* If there have been edits to the main MQLib project since you forked the repository, you should add the main project as an upstream remote to your fork (run `git remote add upstream https://github.com/MQLib/MQLib.git` in the forked repository), and then follow the steps from [here](https://help.github.com/articles/syncing-a-fork/) to sync your fork by fetching the branches from the upstream repository (`git fetch upstream`), checking out your fork's local `master` branch (`git checkout master`), merging with the `upstream/master` branch (`git merge upstream/master`), and pushing the updates to the fork (`git push`).
* On the page for the forked repository on [github.com](https://github.com) (not the main repository), click the "Pull Request" button and then click "Create Pull Request," adding descriptive comments and then clicking "Create Pull Request."
* The package maintainers will review the pull request and add your code and results into the MQLib repository.
