# Instructions to run the code
The files needed to run the code are the load_data.py and cell-count.csv. to run 

'''
python load_data.py
'''

Output files include miraclib_response_boxplots.png and miraclib_response_statistics.csv

Note: Although an AI assistant was used in the completion of this project, I tried to reason through each step and verified the code generated.

# Schema Explanation
The schema was designed to separate different levels of the data rather than storing everything in one large table. This was done to reduce duplication and make the database easier to query.

The core tables were
- projects
- subjects
- samples
- cell types
- cell counts

The projects, subjects, and samples were separate to reduce repitition in the database. For example each patient had 3 samples, but their sex, treatment, condition, and response do not change so those values are stored in subjects.

Integer keys were used because they are compact and efficient for joins and indexes. 

I am not as familiar with SQL and how databases should scale to much larger datasets. From what I can see this structure probably would scale.

## Further analytics
Interesting further comparisions to complete would including lookinng at if there were differences between projects in results, comparing the different cell counts between groups, looking at further effects from age and sex in treatment efficacy. Althrough no significant difference was found in the analysis completed, that was specifically looking at melanoma patients, I did not investigate the other groups. It would also be interesting to see if there were time points that were different between responding and not responding groups. I looked at population counts at time 0, but there might be interesting information in how the cell counts respond to treatment in the other timepoints that might show that the treatment is starting to work or not that would encourage patients to stay on the treatment.

# Code Overview
Right now the code is all in one file. While I know this probably would not be best practice for long term use of this script, my main goal was finishing this submission.

The script has the schema near the top, with analysis functions following, and database generation/loading at the end. 

The format I probably would follow to clean up this code would be to have multiple files with a separation of concerns. The pslit I would follow would be 4 files, load_data.py, database.py, analysis.py, and plots.py. database.py would contain schema creation and query helpers, analysis.py would contain the functions like the analysis section of the current load_data.py, plots.py would contain the code to produce the plots, with load_data.py being the main function and orchestrate the other files. 

## Statistics note
In the statistics I decided against running normality tests so I completed a nonparametric t-test in the Mann-Whitney U test, which compares the two groups of independent samples since I saw no reason the data should be paired. I also added a correction for multiple comparisons in the Benjamin-Hochberg correction since we looked at five populations rather than just one.

# Link to Dashboard
After running
'''
make dashboard
'''

http://localhost:8501