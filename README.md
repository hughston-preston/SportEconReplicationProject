# SportEconReplicationProject
This is a Repository for my Replication Project in Sport Econ 2 

As a part of my Sport Economics II class at Syracuse University, We were tasked to complete a replication project. In this project we were to find a peer-reviewed 
Sport Economics paper published in an academic journal and replicate all that the authors did in their research paper. In addition we were tasked to extend the 
paper's research further by adding an additional element of research that the paper did not include, perhaps another variable of consideration, different years of 
study or new modeling technique. 



For my paper I chose to replicate Kyle Burris and Jacob Coleman's 2018 paper titled "Out of Gas: Quantifying Fatigue in MLB Relievers" that was published Volume 14 
Issue 2 of the Journal of Quantitative Analysis in Sports. This paper seeks to explore the mechanics of fatigue induced velocity effects upon MLB relievers. This is 
accomplished through a statistical approach of Bayesian Inference upon a sequence of pharmacology and toxicology equations that treat pitcher workloads as a fatigue 
inducing toxin that holds diminishing concentrations with each additional day of rest, and possesses unique concentration for each reliever. Therefore, Burris and 
Coleman conclude that certain pitchers are more fatigue resilient while others are more susceptible to velocity decreases as a result of fatigue.

Unfortunately, I am currently not yet knowledgeable of Bayesian statistics and the RJags, RStan, and winBugs packages in R that the Bayesian methodology requires. 
However, Bayesians are merely one interpretation of statistics, and I am currently far better equipped to operate the contrasting frequentist approach towards 
statistics. Therefore, I substitute Burris and Coleman's Bayesian Inference methodology for a frequentist approach based upon linear mixed effects models. I was 
able to maintain a majority of the integrity of Burris and Coleman's initial research process but I do achieve different results than them. I kept their unique 
Pharmacology and toxicology formulas to the greatest extent possible utilizing a frequentist approach.

The additional aspect of research that I completed in my project was a quantile regression that measured how fatigue induced velocity effects differ across the 
reliever velocity spectrum. This is to say I wanted to know whether low velocity pitchers experienced greater or lower fatigue effects than high velocity pitchers.



All of the code was written in R and all of the data was collected using the baseballr package. The R file is listed and entitled Preston_ReplicationProjectCode.R
and the written Replication Paper that covers the project far more thoroughly is listed and entitled PrestonReplicationProject.docx
