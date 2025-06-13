echo "🚀 Starting main pipeline..."

steps=(
  "python main_pipeline_step1.py"
  "python main_pipeline_step2.py"
  "python main_pipeline_step3.py"
  "python main_pipeline_step4.py"
  "python main_pipeline_step5.py"
)

for cmd in "${steps[@]}"; do
  echo -e "\n➡️ Running: $cmd"
  eval "$cmd"
done
