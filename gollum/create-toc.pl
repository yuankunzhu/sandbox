use strict;
use warnings;
use File::Copy;

my $tocBegin = "[//]: # (Start of TOC)\n";
my $tocEnd = "[//]: # (End of TOC)\n";

foreach my $file (<*.md>)
{
    open(my $fh, '<', $file) or die "Can't open $file: $!";
    my @lines = <$fh>;
    close $fh or die "Can't close $file: $!";

    my @headers = ();
    foreach (@lines)
    {
        if ($_ =~ /^###/) {
            push @headers, createLink($_, 3);
        }
        elsif ($_ =~ /^##/) {
            push @headers, createLink($_, 2);
        }
        elsif ($_ =~ /^#/) {
            push @headers, createLink($_, 1);
        }
    }

    if (scalar(@headers) == 0) {
        next;
    }

    open(my $in, '<', $file) or die "Can't open $file' $!";
    open(my $out, '>', "$file.new") or die "Can't write $file.new: $!";

    print $out $tocBegin;
    print $out "**TABLE OF CONTENTS**\n";
    foreach(@headers) {
        print $out $_;
    }
    print $out "\n---\n";
    print $out $tocEnd;

    my $traversingOldToc = 0;
    while(<$in>) {
        if ($_ eq $tocBegin) {
            $traversingOldToc = 1;
        }
        elsif ($_ eq $tocEnd) {
            $traversingOldToc = 0;
        }
        elsif ($traversingOldToc == 0) {
            print $out $_;
        }
    }

    close $in or die "Can't close $file: $!";
    close $out or die "Can't close $file.new: $!";

    move("$file.new", $file) or die "Can't rename $file.new to $file: $!";
}

sub createLink {
    my $currentLine = $_[0];
    my $indent = $_[1];

    my $text = substr($currentLine, $indent);
    $text =~ s/^\s+|\s+$//g;
    my $link = lc $text =~ s/ /-/rg;
    return " " x (($indent-1)*2) . "- " . "<a href=\"#user-content-$link\">$text</a>\n";
}