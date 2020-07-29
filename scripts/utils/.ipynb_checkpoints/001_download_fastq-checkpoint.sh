#! Give your job a name
#SBATCH -J 001_downloading_fastq_files
#! How many cores per task?
#SBATCH --cpus-per-task=1
#! How much memory do you need?
#SBATCH --mem=1G
#! How much wallclock time will be required?
#SBATCH --time=01-00:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL
#! Specify your email address here otherwise you won't recieve emails!
#SBATCH --mail-user=alexander.baker@cruk.cam.ac.uk
#! General partition
#SBATCH -p general

slx_id=$1

mkdir -p data data/fastq

java -jar ../utils/clarity-tools.jar -l $slx_id