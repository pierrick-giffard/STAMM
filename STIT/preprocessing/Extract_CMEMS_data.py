"""
1- You need to have a CMEMS account.
2- In this file, please complete your username and password.
3- Choose your local directory to store data.
4- Choose variables and product to download.
5- Choose the area you want to download and the time range.
6- Choose the prefix of your files.
7- Execute the script.
"""



# Packages
import os
import subprocess
import datetime as dt
import time
import calendar
import sys

try:
    from pathlib import Path
except ImportError:
    print("Trying to Install required module: pathlib\n")
    os.system('python -m pip install pathlib')
try:
    from pathlib import Path
except ImportError:
    print("""You need pathlib module.
                Install it from https://pypi.org/project/pathlib/
                or (recommended) run in a new cell: !pip install pathlib.""")
    
    
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# General Parameters - Tools - Proxy Network - Output Directory
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Path declaration to the motu-client.py opensource-TOOLS to connect to MOTU CopernicusMarineHub.
# If you are not sure, please read the article on "python basic requirements":
# http://marine.copernicus.eu/faq/what-are-the-motu-and-python-requirements/?idpage=169
motu_cl = "python -m motuclient"

# File to log unsuccessful data extraction request(s)
logfile = 'logfile.txt'

# Copernicus Marine API Key - Login Credentials 
username_cmems = 'ldupont'
password_cmems = 'mypassword'

# Proxy Configuration
# Please replace "False" by "True" if you use a proxy to connect to internet and fill in the below variables.
proxy_flag = False
proxy_server_url = "http://your_proxy_url.com"
proxy_server_port = "port"
proxy_user_login = "your_proxy_user_login"
proxy_user_password = "your_proxy_user_password"

# Output directory name to store/save the Copernicus Marine data
local_storage_directory_name = '/data/PSY4/'

# - - - - - - - - - - - - - - - - - - - - - - - - -
# Product(s), Dataset(s) and MOTU server Parameters 
# - - - - - - - - - - - - - - - - - - - - - - - - -

''' /!\ All Copernicus marine datasets are NOT hosted by a single server
        You can always rely on the "VIEW SCRIPT" button of the Copernicus Marine Website (marine.copernicus.eu),
        using its DataExtraction WebInterface (also called GUI).
        It will generate the parameters of your extraction settings based on your selection.
        Please refer to this article to understand how to call/trigger this webservice/feature to generate the right parameters:
        http://marine.copernicus.eu/faq/how-to-write-and-run-the-script-to-download-cmems-products-through-subset-or-direct-download-mechanisms/?idpage=169'''


# set ocean variables in a list separated by a comma, as follow:
# variables = ['variable1', 'variable2', etc]
variables = ['thetao', 'uo', 'vo']

# set the product identifier name, its dataset identifier name, a data request identifier name (can be anything) and the motu server URL
product_name       = 'GLOBAL_ANALYSIS_FORECAST_PHY_001_024'
dataset_name       = 'global-analysis-forecast-phy-001-024'
data_request_name  = 't-u-v-dailymean'
motu_server        = 'http://nrtcmems.mercator-ocean.fr/motu-web/Motu'


# - - - - - - - - - - - - - - - - - - - - - -
# Geographical Area Parameters and Timerange
# - - - - - - - - - - - - - - - - - - - - - -

#  -y LATITUDE_MIN, --latitude-min=LATITUDE_MIN
#                        The min latitude (float in the interval [-90 ; 90])
#  -Y LATITUDE_MAX, --latitude-max=LATITUDE_MAX
#                        The max latitude (float in the interval [-90 ; 90])
#  -x LONGITUDE_MIN, --longitude-min=LONGITUDE_MIN
#                        The min longitude (float in the interval [-180 ; 180])
#  -X LONGITUDE_MAX, --longitude-max=LONGITUDE_MAX
#                        The max longitude (float in the interval [-180 ; 180])
#  -z DEPTH_MIN, --depth-min=DEPTH_MIN
#                        The min depth (float in the interval [0 ; 2e31] or
#                        string 'Surface')
#  -Z DEPTH_MAX, --depth-max=DEPTH_MAX
#                        The max depth (float in the interval [0 ; 2e31] or
#                        string 'Surface')


# Geographical Area (x west-east longitude, y south-north latitude, z depth)
# set the bounding box(es) of interest, as follow: 
#     boundingboxID = [minLon, maxLon, minLat, maxLat, minSurface, maxSurface]
# /!\ the order of parameters is important (starting with minimum longitude and finishing with maximum depth)
boundingbox1  = [-180,179.91,-80,90,0.49,0.50]

# Date - Timerange
yyyystart = 2009
mmstart = 1
yyyyend = 2018
mmend = 12
hhstart = " 12:00:00"
hhend = " 12:00:00"
dd = 31
DaylySerie = True #option to write daily files.

# Output prefix file name
pre_name = "CMEMS_PSY4_"


# Check motuclient 
try:
    import motuclient
except ImportError:
    print('\nThe module motuclient was not found. To overcome this problem, please:\n- Open a terminal and type in: pip install motuclient\nOr\n- Create a new cell and write in: !pip install motuclient')
    sys.exit()

# Check output directory, create if necessary
p = Path(local_storage_directory_name)
if not p.exists():
    p.mkdir(parents=True)

# Flags to let the server clears the buffer - better to be respectful when retrieving OPEN and FREE data
buffer_flag = False
cmd_flag = False

# Error Handle on dates (to illustrate an if statement with eval param '>')
if yyyystart > yyyyend:
    print("[ERROR] in [Date Parameters]")
    print("""Please double check your date parameters, specifically the "yyyystart" which is currently greater than "yyyyend.""")
    print("""End of data extraction service.""")
    sys.exit()

# Other variable definitions to be compatible with deprecated script versions still available on the Internet
local_storage_directory_name = Path(local_storage_directory_name)
logfilepath = local_storage_directory_name.absolute() / logfile
log_cmems = "-u " + username_cmems
pwd_cmems = "-p " + password_cmems
pre_fic_cmd = "-f " + pre_name
out_cmd = "-o " + str(local_storage_directory_name.absolute())
proxy_user = "--proxy-user " + proxy_user_login
proxy_pwd = "--proxy-pwd " + proxy_user_password
proxy_server = "--proxy-server " + proxy_server_url + ":" + proxy_server_port

# To illustrate a simple Error Handle to delete a file when desired
try:
    Path.unlink(logfilepath)
except OSError:
    print ("")

# Define a dict based on the above boundingboxes for the geographical regions of interest of the `data_request_name`
#dict_boundingbox = {'boundingbox1':boundingbox1, 'boundingbox2':boundingbox2, 'boundingbox3':boundingbox3, 'boundingbox4':boundingbox4, 'boundingbox5':boundingbox5}
dict_boundingbox = {'boundingbox1':boundingbox1}


variable_string=' --variable '
all_variables = ''
for vv in variables:
    all_variables = variable_string + vv + all_variables
    
dict_dataset = {data_request_name: [all_variables, "--product-id " + dataset_name, "--service-id " + product_name + "-TDS", "--motu " + motu_server]}

print("\n+----------------------------+\n| ! - CONNEXION TO CMEMS HUB |\n+----------------------------+\n\n")
    
# To illustrate a For_Loop in order to generate download requests for several datasets held in a product
for k, v in dict_dataset.items():
    for key, value in dict_boundingbox.items():
        minX = value[0]
        maxX = value[1]
        minY = value[2]
        maxY = value[3]
        minZ = value[4]
        maxZ = value[5]
        list_abs = [abs(maxX-minX),abs(maxY-minY),abs(maxZ-minZ)]
        var     = v[0]
        dataset = v[1]
        service = v[2]
        motu    = v[3]
        
        if buffer_flag:
            print ("Little pause to let the server clearing the buffer, it will AUTOMATICALLY resume once it's completed.\nNot mandatory but server-friendly :-)\n")
            time.sleep(5)
            buffer_flag = False

        # Date declaration
        date_start = dt.datetime(yyyystart,mmstart,dd,0,0)
        date_end = dt.datetime(yyyyend,mmend,dd,0,0)

        # To illustrate a While_Loop in order to extract dailymean data, packed by year, 
        # for as many download requests as number of years available in the timerange.
        while (date_start <= date_end):
            # To illustrate a check on applying a stack (either by month or year) based on the number of gridpoints.
            # It's known that there is a higher chance of getting an error if too many gridpoints for stack by year.
            if DaylySerie:
                date_end_cmd = date_start
                ficout = pre_name + k + "_daystack_" + key + "_" + date_end_cmd.strftime("%Y-%m-%d") + ".nc"
                
            elif len([i for i in list_abs if i > 40]) > 0:
                date_end_cmd = (dt.datetime(date_start.year, date_start.month,\
                                            calendar.monthrange(date_start.year, date_start.month)[1]))
                ficout = pre_name + k + "_monthstack_" + key + "_" + date_end_cmd.strftime("%Y-%m") + ".nc"
                
            else:
                date_end_cmd = (dt.datetime(date_start.year, date_end.month,\
                                            calendar.monthrange(date_start.year, date_end.month)[1]))
                ficout = pre_name + k + "_yearstack_" + key + "_" + date_end_cmd.strftime("%Y") + ".nc"
            date_cmd = ' -t \"' + date_start.strftime("%Y-%m-%d") + hhstart + '\"' + ' -T \"' + date_end_cmd.strftime("%Y-%m-%d") + hhend + '\"'
            fic_cmd = '-f ' + ficout
            print("----------------------------------\n- ! - Processing dataset request : %s\n----------------------------------\n"%ficout)
            if not Path(local_storage_directory_name / ficout).exists():
                if proxy_flag:
                    if minZ is None:
                        cmd = ' '.join([motu_cl, log_cmems, pwd_cmems,\
                                    motu, service, dataset,\
                                    "-x", str(minX), "-X", str(maxX), "-y", str(minY), "-Y", str(maxY),\
                                    date_cmd, var, out_cmd, fic_cmd,\
                                    proxy_server, proxy_user, proxy_pwd, "-q"])
                    else:
                        cmd = ' '.join([motu_cl, log_cmems, pwd_cmems,\
                                    motu, service, dataset,\
                                    "-x", str(minX), "-X", str(maxX), "-y", str(minY), "-Y", str(maxY), "-z",str(minZ), "-Z", str(maxZ),\
                                    date_cmd, var, out_cmd, fic_cmd,\
                                    proxy_server, proxy_user, proxy_pwd, "-q"])
                else:
                    if minZ is None:
                        cmd = ' '.join([motu_cl, log_cmems, pwd_cmems,\
                                    motu, service, dataset,\
                                    "-x", str(minX), "-X", str(maxX), "-y", str(minY), "-Y", str(maxY),\
                                    date_cmd, var, out_cmd, fic_cmd, "-q"])
                    else:
                        cmd = ' '.join([motu_cl, log_cmems, pwd_cmems,\
                                    motu, service, dataset,\
                                    "-x", str(minX), "-X", str(maxX), "-y", str(minY), "-Y", str(maxY), "-z",str(minZ), "-Z", str(maxZ),\
                                    date_cmd, var, out_cmd, fic_cmd, "-q"])
                print("## MOTU API COMMAND ##")
                print(cmd)
                print("\n[INFO] CMEMS server is checking both your credentials and command syntax. If successful, it will extract the data and create your dataset on the fly. Please wait. \n")
                subpro=subprocess.Popen(cmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                message,erreur = subpro.communicate()
                stat = subpro.returncode
                if stat != 0:
                        print("-- ERROR Incorrect Credentials :\n %s"%message)
                        with open(logfilepath,'a') as mylog:
                            mylog.write("Error : %s NOK\nDue to : %s"%(ficout,message))
                        if b'HTTP Error 400' in message:
                            print('HTTP Error 400')
                            sys.exit()
                        if b'HTTP Error 407' in message:
                            print ('''[INFO] Proxy Authentication Required to connect to the Central Authentication System https://cmems-cas.cls.fr/cas/login\n\n[INFO] Check the value of proxy_flag (it should be True).\n\n[INFO] Double check your proxy settings:\n  --proxy-server=PROXY_SERVER\n                        the proxy server (url)\n  --proxy-user=PROXY_USER\n                        the proxy user (string)\n  --proxy-pwd=PROXY_PWD\n                        the proxy password (string)\n\n[INFO] If your proxy credentials are correct but your proxy password (string) contains a '@' then replace it by '%%40' ''')
                            print ('''[INFO] This issue is raised due either a misconfiguration in proxy settings or a network issue. If it persists, please contact your network administrator.''')
                            sys.exit()
                        print("""[INFO] Failed data extraction has been logged.\n""")
                else:
                    if b"[ERROR]" in message:
                        print("-- ERROR Downloading command :\n %s"%message)
                        with open(logfilepath,'a') as mylog:
                            mylog.write("Error : %s NOK\nDue to : %s"%(ficout,message))
                        print ("""[INFO] Failed data extraction has been logged.\n""")
                    else:
                            print("-- MOTU Download successful :\n %s OK\n"%fic_cmd.split()[1])
                            print("5 seconds break to let the server releasing the token to grant next request. It will resume AUTOMATICALLY.\n")
                            time.sleep(5)
                            cmd_flag = True
            else:
                print("-- Your dataset for %s has already been downloaded in %s --\n"% (fic_cmd.split()[1],out_cmd.split()[1]))
                cmd_flag = False

            date_start = date_end_cmd + dt.timedelta(days=1)

        if cmd_flag:
            buffer_flag = True
            cmd_flag = False
    
if not logfilepath.exists():
    print("\n------------------------------------------------\n - ! - Your Copernicus Dataset(s) are located in %s\n------------------------------------------------\n"%(out_cmd.split()[1]))
else :
    print("## [ERROR] ##")
    print("/!\\ Some download requests failed. Please see recommendation in %s."%(logfilepath))
    print("/!\\ Do not move netCDF files and relaunch the script to automatically download only failed data request(s) and to skip successful.")
print("+--------------------------------------------+\n| ! - CONNEXION TO CMEMS HUB HAS BEEN CLOSED |\n+--------------------------------------------+\n")
#------------------------------------------------- End of Script -----------------------------------------------------