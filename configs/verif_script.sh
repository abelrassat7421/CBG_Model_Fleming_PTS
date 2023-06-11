for ((i=601; i<=1200; i++)); do
    if [ ! -f "pts_test_${i}.yml" ]; then
        echo "File pts_test_${i}.yml does not exist."
    fi
done

