/*----------------------------------------------------------------------
 2005.07.25 sample_decoder.c

 GRIB2 ファイル サンプルデコード処理プログラム

 このプログラムの全部又は一部を利用してもかまいませんが、利用したこ
 とによって利用者が被った直接的又は間接的ないかなる損害についても、
 気象庁は一切責任を負いません。
 また、プログラムに関する個別の対応は行いかねますので、ご容赦願います。

 利用方法

 ANSI 準拠の C コンパイラでコンパイルして下さい。
 GRIB2 ファイルのファイル名を引数に与えて実行すると、
 ファイルの内容が表示されます。
 ----------------------------------------------------------------------*/

 # include <stdio.h>
 # include <stdlib.h>
 # include <math.h>


 /* 1,2,4,8 : signed integer (x byte)
          u : unsigned integer (1 byte)
          S : unsigned integer (2 byte)
          C : character (4 byte)
          R : float (4 byte) */

 char *sectionFormat[9] = {
    "C2uu8", /* Section 0 */
    "4u22uuu2uuuuuuu", /* Section 1 */
    "", /* Section 2 */
    "4uu4uuS", /* + Section 3 */
    "4u2S", /* + Section 4 */
    "4u4S", /* + Section 5 */
    "4uu", /* + Section 6 */
    "4u", /* + Section 7 */
    "C" /* Section 8 */
 };

 typedef struct { int secno, templat; char *format; } TemplateFormat;

 TemplateFormat templateFormat[] = {
    { 3, 0, "uu4u4u4444444u4444u" },
    { 4, 0, "uuuuu2uu4u14u14" },
    { 4, 1, "uuuuu2uu4u14u14uuu" },
    { 4, 8, "uuuuu2uu4u14u142uuuuuu4uuu4u4" },
    { 4,11, "uuuuu2uu4u14u14uuu2uuuuuu4uuu4u4" },
    { 5, 0, "R22uu" },
    { 0, 0, "" }
 };

 static int isLittleEndian;

 # define FIT_BYTE_ORDER(pointer,size) if(isLittleEndian)swabN(pointer,size,1)

 /* reverse byte order */
 static void
 swabN( void *buf, int size, int nn )
 {
    char *ba, *bb, *buf2 = buf; 
    while( nn-- ) {
        bb = ( ba = buf2 ) + size -1;
        do {
            char a;
            a = *ba;
            *ba = *bb;
            *bb = a;
        } while( ++ba < --bb );
      buf2 += size;
    }
 }

 int
 read_section_0( FILE *fp, void **sec_buffer )
 {
    char *bufr;
    *sec_buffer = realloc( *sec_buffer, 16 );
    bufr = *sec_buffer;
    if( fread( bufr, 1, 16, fp ) != 16 ||
    strncmp( bufr, "GRIB", 4 ) != 0 ) {
        fprintf( stderr, "Really GRIB file ?¥n" );
        exit(1);
    }
    return 0;
 }

 int
 read_section_X( FILE *fp, void **sec_buffer )
 {
    int length;
    unsigned char *bufp;

    fread( &length, 4, 1, fp );
    if( strncmp( (char *)(&length), "7777", 4 ) == 0 ) {
        *sec_buffer = realloc( *sec_buffer, 4 );
        strncpy( *sec_buffer, "7777", 4 );
        return 8;
    }
    FIT_BYTE_ORDER( &length, 4 );

    *sec_buffer = realloc( *sec_buffer, length );
    bufp = *sec_buffer;
    if( fread( bufp + 4, 1, length - 4, fp ) != length - 4 ) {
        fprintf( stderr, "Unexpected EOF¥n" );
        exit(1);
    }
    FIT_BYTE_ORDER( &length, 4 );
    memcpy( bufp, &length, 4 );

    return (int)bufp[4]; /* section No. */
 }

 void
 decode_buf( char *sec_buffer, int *index, char *format, double **dvpp )
 {
    unsigned char buffer[8];
    int size, ii, jj, missing;
    double *dvp, dd;

    dvp = *dvpp;
    
    for( ii = *index ; *format ; format++ ) {
        switch( *format ) {
            case '1' :
            case '2' :
            case '4' :
            case '8' : size = *format - '0'; break;
            case 'u' : size = 1; break;
            case 'S' : size = 2; break;
            case 'R' :
            case 'C' : size = 4; break;
            default :
            fprintf( stderr, "Internal Error !¥n" );
            exit(9);
        }
        memcpy( buffer, sec_buffer+ii, size );

        missing = 1;
        for( jj = 0 ; jj < size ; jj++ )
        missing &= ( buffer[jj] == 0xFF );

        switch( *format ) {
            case 'u' : dd = *buffer; break;
            case 'S' : dd = (buffer[0]<<8) + buffer[1]; break;

            case '1' :
            case '2' :
            case '4' :
            case '8' : dd = buffer[0] & 0x7F;
            for( jj = 1 ; jj < size ; jj++ )
            dd = dd * 256.0 + buffer[jj];
            if( *buffer & 0x80 ) dd = -dd; break;

            case 'R' : FIT_BYTE_ORDER( buffer, size );
            dd = *(float *)buffer; break;

            case 'C' : memcpy( &dd, buffer, size ); break;
        }
        if( size == 1 ) printf( " %4d : ", ii+1 );
        else printf( "%4d .. %-4d: ", ii+1, ii+size );

        if( missing ) printf( "missing ¥n" );
        else if( *format == 'R' ) printf( "%f¥n", dd );
        else if( *format == 'C' ) printf( "'%.4s'¥n", &dd );
        else printf( "%.0f¥n", dd );

        *dvp++ = dd;
        ii += size;
    }
    *index = ii;
    *dvpp = dvp;
}

void
decode_section( int secno, char *sec_buffer, double **double_values )
{
    int index, length;
    double *dvp;
    char *format;
    TemplateFormat *tfp;

    switch( secno ) {
    case 0 : length = 5; break;
    case 1 : length = 15; break;
    case 7 : length = 2; break;
    case 8 : length = 1; break;
    default :
    memcpy( &length, sec_buffer, 4 );
    FIT_BYTE_ORDER( &length, 4 );
    }
    dvp = realloc( *double_values, sizeof(double) * length );
    *double_values = dvp;

    printf( "== SECTION %d == ¥n", secno );

    format = sectionFormat[secno];
    index = 0;

    decode_buf( sec_buffer, &index, format, &dvp );

    if( 3 <= secno && secno <= 5 ) {
        int templat = (int)dvp[-1];
        for( tfp = templateFormat ; tfp->secno ; tfp++ ) {
            if( tfp->secno == secno &&
             tfp->templat == templat ) break;
        }
        if( ! tfp->secno ) {
            fprintf( stderr, "No Information about Template %d.%d¥n",
            secno, templat );
            return;
        }
        decode_buf( sec_buffer, &index, tfp->format, &dvp );
    }
 }

 void
 unpack_data( const char *sec_7, const double *values_5, float **out )
 {
    double rr, ee, dd, pow_2_e, pow_10_d;
    float *op;
    int num, nbit, ii;
    unsigned int uw, mask;
    const unsigned char *uc;

    uc = (const unsigned char *)sec_7 + 5;

    num = values_5[2];
    rr = values_5[4];
    ee = values_5[5];
    dd = values_5[6];
    nbit = values_5[7];
    pow_2_e = pow( 2.0, ee);
    pow_10_d = pow(10.0, dd);

    *out = realloc( *out, sizeof(**out) * num );
    op = *out;

    mask = 1;
    for( ii = 0 ; ii < nbit - 1 ; ii++ )
    mask = mask << 1 | 1;

    for( ii = 0 ; ii < num ; ii++ ) {
        memcpy( &uw, &uc[(nbit*ii)/8], 4 ); 
        FIT_BYTE_ORDER( &uw, 4 );
        uw = ( uw>>(32-nbit-(nbit*ii)%8)) & mask;
        *op++ = ( rr + pow_2_e * uw ) / pow_10_d;
    }
 }

 # include <float.h>

 char *
 elem_name( const double *sec4 )
 {
    static char name_bufr[256];
    char level[256], elem[256], *elemp;
    int cc, nn, tt;

    typedef struct { int category, number; char *name; } ETable;

    static ETable etable [] = {
        { 0, 0, "Temperature" },
        { 1, 1, "RelativeHumidity" },
        { 1, 8, "TotalPrecipitation" },
        { 2, 2, "UofWind" },
        { 2, 3, "VofWind" },
        { 2, 8, "VerticalVelocity(Pressure)" },
        { 3, 0, "Pressure" },
        { 3, 1, "PressureReducedToMSL" },
        { 3, 5, "GeopotentialHeight" },
        { 6, 1, "TotalCloudCover" },
        { 6, 3, "LowCloudCover" },
        { 6, 4, "MediumCloudCover" },
        { 6, 5, "HighCloudCover" },
        { 0, 0, NULL }
    };
    ETable *etp;

    cc = sec4[4];
    nn = sec4[5];
    for( etp = etable ; etp->name ; etp++ )
    if( etp->category == cc && etp->number == nn ) {
        elemp = etp->name;
        break;
    }
    if( ! etp->name ) {
        fprintf( stderr, "No information about Product category %d number %d¥n",
        cc, nn );
        sprintf( elem, "Cat(%d)num(%d)", cc, nn );
        elemp = elem;
    }
    if( sec4[13] == 1 ) strcpy( level, "Surface" );
    else if( sec4[13] == 101 ) strcpy( level, "MeanSeaLevel" );
    else if( sec4[13] == 103 )
        sprintf( level, "%.1fm", sec4[15] / pow(10., sec4[14]) );
    else if( sec4[13] == 100 )
        sprintf( level, "%.0fhPa", sec4[15] / pow(10.,sec4[14])/100. );
    else sprintf( level, "UnknownFixedsurface(%.0f)", sec4[13] );

    switch( (int)sec4[3] ) {
        case 8 : tt = sec4[12] + sec4[30]; break;
        case 11 : tt = sec4[12] + sec4[33]; break;
        default : tt = sec4[12]; break;
    } 
    switch( (int)sec4[3] ) {
        char *ensm;
        case 1: case 11:
        switch( (int)sec4[19] ) {
            case 0 : ensm = "H"; break;
            case 1 : ensm = "" ; break;
            case 2 : ensm = "m"; break;
            case 3 : ensm = "p"; break;
        }
        sprintf( name_bufr, "%s_%s_T%d_M%.0f%s",
        elemp, level, tt, sec4[20], ensm ); break;
        default:
        sprintf( name_bufr, "%s_%s_T%d", elemp, level, tt ); break;
    }
    return name_bufr;
 }

 void
 show_data_statistics( const double *values_4, const float *unpacked_data, int num )
 {
    int ii;
    double sum = 0, min = DBL_MAX, max = -DBL_MAX;

    printf( " << %s >> ¥n", elem_name( values_4 ) );

    for( ii = 0 ; ii < num ; ii++ ) {
        double ff = unpacked_data[ii];
        sum += ff;
        if( min > ff ) min = ff;
        if( max < ff ) max = ff;
    }
    printf( " ( num = %d, min = %f, max = %f, average = %f )¥n",
            num, min, max, sum / num );
 }

 void
 save_float_file( const double *values_4, const float *unpacked_data, int num )
 {
    FILE *fpout;
    char fileName[256];

    sprintf( fileName, "%s.dat", elem_name( values_4 ) );

    fprintf( stderr, "output '%s'¥n", fileName );

    if( ( fpout = fopen( fileName, "wb" ) ) == NULL ) {
        fprintf( stderr, "file '%s' open error! ¥n", fileName );
        exit(1);
    }
    fwrite( unpacked_data, sizeof(float), num, fpout );
    fclose( fpout );
 }

 int main( int argc, char *argv[] )
 {
    FILE *fpin;
    char *fileName;
    int secno, ii;
    void *sec_buffer = NULL;
    float *unpacked_data = NULL;
    double *sec_double_value[9]; 


    /* check system byte order */
    isLittleEndian = 1;
    isLittleEndian = *(char *)(&isLittleEndian);

    for( ii = 0 ; ii < sizeof(sec_double_value)/sizeof(*sec_double_value) ; )
    sec_double_value[ii++] = NULL;
    if( argc == 1 ) {
        fprintf( stderr, "¥n¥n usage: %s 'grib2 file name'¥n¥n", argv[0] );
        exit(1);
    }
    fileName = argv[1];
    if( ( fpin = fopen( fileName, "rb" ) ) == NULL ) {
        fprintf( stderr, "grib2 file '%s' open error! ¥n", fileName );
        exit(1);
    }
    secno = read_section_0( fpin, &sec_buffer );
    decode_section( secno, sec_buffer, &sec_double_value[secno] );

    while( secno = read_section_X( fpin, &sec_buffer ) ) {
        decode_section( secno, sec_buffer, &sec_double_value[secno] );
        if( secno == 8 ) break;
        if( secno == 7 ) {
            unpack_data( sec_buffer, sec_double_value[5], &unpacked_data );

            show_data_statistics( sec_double_value[4], unpacked_data,
            (int)sec_double_value[5][2] );
            /*
            save_float_file ( sec_double_value[4], unpacked_data,
            (int)sec_double_value[5][2] );
            */
        }
    }
    fclose(fpin);
    return 0;
 } 