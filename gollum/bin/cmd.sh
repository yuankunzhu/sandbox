git clone git@github.com:your-account/your-project.wiki.git
brew install icu4c
sudo gem install gollum # willl hang on "Installing ri documentation for gollum-4.0.1" for a while
git init & gollum

# covert md to html
brew install pandoc
pandoc -f markdown some-file.md > some-file.html
