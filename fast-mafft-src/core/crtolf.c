/* fast-mafft (overhead lever): the mafft script copies its input with
       cat "$1" | tr "\r" "\n" > "$TMPFILE/infile"
   and tr only copies about 30-50 MB/s. For a file whose bytes are all ASCII (< 0x80), tr gives
   the same bytes in every locale (C, UTF-8, SJIS, eucJP): each CR becomes LF and everything
   else, NUL and other control bytes included, is copied unchanged. This program does that
   copy. If the file has any byte >= 0x80, where tr's result depends on the locale (a UTF-8
   locale rejects invalid sequences, for example), it stops with exit status 3 and the script
   runs the original tr pipeline instead, which rewrites the output file.

   usage: crtolf FILE > OUT    (exit 0: done; 2: I/O error; 3: non-ASCII byte found) */
#include <stdio.h>

int main( int argc, char *argv[] )
{
	static unsigned char buf[1<<20];
	size_t n, i;
	FILE *fp;

	if( argc != 2 ) return( 2 );
	fp = fopen( argv[1], "rb" );
	if( fp == NULL ) return( 2 );
	while( ( n = fread( buf, 1, sizeof( buf ), fp ) ) > 0 )
	{
		for( i=0; i<n; i++ )
		{
			if( buf[i] >= 0x80 ) return( 3 );
			if( buf[i] == '\r' ) buf[i] = '\n';
		}
		if( fwrite( buf, 1, n, stdout ) != n ) return( 2 );
	}
	if( ferror( fp ) ) return( 2 );
	fclose( fp );
	if( fflush( stdout ) != 0 || ferror( stdout ) ) return( 2 );
	return( 0 );
}
