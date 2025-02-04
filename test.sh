for i in {30..1}; do
    GIT_AUTHOR_DATE="$(date --date="$i days ago" "+%Y-%m-%dT12:00:00")"
    GIT_COMMITTER_DATE="$GIT_AUTHOR_DATE"
    echo "Update on day $i" > fake_update.txt
    git add .
    git commit -m "Regular update on day $i" --date="$GIT_AUTHOR_DATE"
done
git push origin main
