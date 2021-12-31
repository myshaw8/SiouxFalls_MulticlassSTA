This notebook requires additional dependencies that are not initially included in network-modeling.
In order to update the anaconda notebook environment, please follow these steps:

1. Open Anaconda Prompt.
2. Type "conda env update -n network-modeling -f " (no quotations).
3. Find the full path to the "network-modeling_environment.yml" file inside the project folder.
4. Add the path to the end of the command. It should include the environment config file.
5. Your full command should look something like this:

conda env update -n network-modeling -f C:\Users\MYShaw\Documents\Jupyter\ISYE8803\Team Homework\network-modeling_environment.yml

6. Press enter to execute the command. Your packages should update. You may need to hit 'y' to accept the installation of additional packages.