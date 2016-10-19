#!/bin/bash
# Based on
# http://sleepycoders.blogspot.com.au/2013/03/sharing-travis-ci-generated-files.html
# https://help.github.com/articles/creating-project-pages-manually
#
# This differs by deleting the history of the gh pages branch by doing
# a force push off a brand new branch.  This helps keep the size of
# the repo under control.
echo -e "Preparing to copy generated files to gh-pages branch"
#echo $TRAVIS_PULL_REQUEST
#echo $TRAVIS_BRANCH
#echo $TRAVIS_R_VERSION
if [[ "$TRAVIS_PULL_REQUEST" == 'false' && "$TRAVIS_BRANCH" == 'master' ]]; then
        echo -e "Starting to update gh-pages\n"

        mkdir -p $HOME/docs

        if $PKGDOWN_BUILT; then

          #go to home and setup git
          cd $HOME
          git config --global user.email "ross@ecohealthalliance.org"
          git config --global user.name "Noam Ross (bot)"

          echo -e "Recloning project"
          # Reclone the project, using the secret token.  Uses /dev/null to avoid leaking decrypted key
          git clone --quiet --branch=gh-pages --single-branch https://${GH_TOKEN}@github.com/ecohealthalliance/metaflu.git gh-pages > /dev/null

          cd gh-pages

          # Move the old branch out of the way and create a new one:
          git branch -m gh-pages-old
          git checkout --orphan gh-pages

          # Delete all the files and replace with our good set
          git rm -rf .
          cp -Rf $HOME/docs/* .

          # add, commit and push files
          git add -f .
          git commit -m "Travis build $TRAVIS_BUILD_NUMBER pushed to gh-pages"
          echo -e "Pushing to origin/gh-pages"
          git push -fq origin gh-pages > /dev/null

          echo -e "Uploaded generated files to gh-pages\n"
        else
          echo -e "Document building failed; not uploading website."
        fi
else
    echo -e "Not pushing pages - only do this for master branch updates"
fi
