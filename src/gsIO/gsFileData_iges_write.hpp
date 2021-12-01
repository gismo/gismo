
char o_crea_ent = 's';        // p, g, c, s for corresponding entity to write
int  entity_coords[MAXENTITY]; // A number of coordinates in entity */
int  o_incr_layer = 0;        // 1 = increase layer */

void fill_int(int n, int poz, int x, char *s2)
{
    /* fill "x"-characters array from position "poz" in string "s2" with the value of "n" (right aligned) */
    int i, j;
    char s1[81];
    int s1_len;
    sprintf(s1, "%d", n);
    s1_len=strlen(s1);
    j=poz;
    while (j < (poz+x-s1_len)) {
        s2[j]=' ';
        j++;
    }
    i=0;
    j=poz+x-s1_len;
    while (i < s1_len) {
        s2[j]=s1[i];
        i++;
        j++;
    }
}

void fill_s(const char *s1, int poz, int x, char *s2)
{
    /* fill "x"-characters array from position "poz" in string "s2" with the string "s1" (right aligned) */
    int i, j;
    int s1_len;
    s1_len=strlen(s1);
    j=poz;
    while (j < (poz+x-s1_len)) {
        s2[j]=' ';
        j++;
    }
    i=0;
    j=poz+x-s1_len;
    if (s1_len > x) s1_len=x; /* just to be sure !*/
    while (i < s1_len) {
        s2[j]=s1[i];
        i++;
        j++;
    }
}

void writeiges_de(int ** entity, int & emark, int & entity_sum, FILE * fw)
{
    /* write DE section to IGES file */
    char first_row[81], second_row[81];
    /*[0     ][8     ][16    ][24    ][32    ][40    ][48    ][56    ][64    ][72    ]*/
    strcpy(first_row,   "     xxx       x       0       0       x       0       0       000000000D      x");
    strcpy(second_row,  "     xxx       0       1       x       x       0       0       x       0D      x");
    char const * cc = "3d point";
    switch (o_crea_ent)
    {
    case 'p':
        fill_int(116, 0, 8, first_row);                  /* 116 - point */
        fill_int(116, 0, 8, second_row);
        fill_int(0, 32, 8, second_row);
        for (emark=0; emark < entity_sum; emark++) {
            fill_int((2*emark+1), 73, 7, first_row);
            fill_int((2*emark+2), 73, 7, second_row);
            cc = "3d point";
            fill_s(cc, 56, 8, second_row);
            fill_int(entity[emark][PD_PTR], 8, 8, first_row);
            fill_int(entity[emark][PD_CNT], 24, 8, second_row);
            fill_int(entity[emark][ELAYER], 32, 8, first_row);
            fprintf(fw, "%s\n", first_row);
            fprintf(fw, "%s\n", second_row);
        }
        break;
    case 'g':
        fill_int(106, 0, 8, first_row);                  /* 106 - group of points */
        fill_int(106, 0, 8, second_row);
        fill_int(2, 32, 8, second_row);
        for (emark=0; emark < entity_sum; emark++) {
            fill_int((2*emark+1), 73, 7, first_row);
            fill_int((2*emark+2), 73, 7, second_row);
            fill_s("pt_group", 56, 8, second_row);
            fill_int(entity[emark][PD_PTR], 8, 8, first_row);
            fill_int(entity[emark][PD_CNT], 24, 8, second_row);
            fill_int(entity[emark][ELAYER], 32, 8, first_row);
            fprintf(fw, "%s\n", first_row);
            fprintf(fw, "%s\n", second_row);
        }
        break;
    case 'c':
        fill_int(126, 0, 8, first_row);                  /* 126 - B-spline curve */
        fill_int(126, 0, 8, second_row);                 /* 126 */
        fill_int(0, 32, 8, second_row);                  /* form (subtype) */
        for (emark=0; emark < entity_sum; emark++) {
            fill_int((2*emark+1), 73, 7, first_row);       /* DE seq nr */
            fill_int((2*emark+2), 73, 7, second_row);      /* DE seq nr */
            fill_s("3d BsCrv", 56, 8, second_row);         /* entity label */
            fill_int(entity[emark][PD_PTR], 8, 8, first_row);   /* PD pointer */
            fill_int(entity[emark][PD_CNT], 24, 8, second_row); /* PD count */
            fill_int(entity[emark][ELAYER], 32, 8, first_row);  /* layer */
            fprintf(fw, "%s\n", first_row);
            fprintf(fw, "%s\n", second_row);
        }
        break;
    case 's':
        fill_int(128, 0, 8, first_row);                  /* 128 - B-spline curve */
        fill_int(128, 0, 8, second_row);                 /* 128 */
        fill_int(0, 32, 8, second_row);                  /* form (subtype) */
        for (emark=0; emark < entity_sum; emark++) {
            fill_int((2*emark+1), 73, 7, first_row);       /* DE seq nr */
            fill_int((2*emark+2), 73, 7, second_row);      /* DE seq nr */
            fill_s("3d BsSrf", 56, 8, second_row);         /* entity label */
            fill_int(entity[emark][PD_PTR], 8, 8, first_row);   /* PD pointer */
            fill_int(entity[emark][PD_CNT], 24, 8, second_row); /* PD count */
            fill_int(entity[emark][ELAYER], 32, 8, first_row);  /* layer */
            fprintf(fw, "%s\n", first_row);
            fprintf(fw, "%s\n", second_row);
        }
        break;
    }
}

void pd_linestring_complement(char *s, int ptr_de, int seq) {
    int i;
    char endstring[20];
    for (i=strlen(s); i <= 64; i++) {
        s[i]=' ';
    }
    s[65]='\0';
    sprintf(endstring, "%07dP%7d", ptr_de, seq);
    strcat(s, endstring);
}

int is_numberlike(char znak) {
    if ( (znak >= '0' && znak <= '9')
         || znak == '-'
         || znak == '+'
         || znak == '.'
         || znak == ','
         || znak == 'E'
         || znak == 'e'
         || znak == 'D'
         || znak == 'd' ) {
        return 1;
    }
    else {
        return 0;
    }
}

void writeiges_pd(int ** entity, int & emark, int & entity_sum, FILE * fw)
{
    /* write PD section to IGES file */
    int i;
    char znak;
    int  i_znak;

    char riadok[81];
    char zbytok[81];

    char x[FIELD_L+1], y[FIELD_L+1], z[FIELD_L+1];
    char *xyz;      /* pointer to x[], y[] or z[] */
    int seq;        /* PD line position (seq number) */
    int line_count; /* number of coordinates for current entity */

    /* State variable for 116 and 106 */
    int complet;

    /* State variables for 126 */
    int prev_char_numberlike;               /* 0,1  Is the previous character number-like? */
    int curr_char_numberlike;               /* 0,1  Is the current character number-like? */
    int numberlike_in_line_counter;         /* 0,1,2,... How many characters in line are number-like? */
    int the_last_coord; /* last coord being started (0= no coordinate yet, xyz=123) */
    int prev_line_nonblank;                 /* 0,1  Contains the previous line number-like characters? */
    int curr_line_nonblank;                 /* 0,1  Contains the current line number-like characters? */

    /* Initializers - global */
    emark= 0;

    FILE * fr=NULL;//REMOVE (input txt file)
    fseek(fr, 0L, SEEK_SET);

    /* Initializers - local */
    znak= ' ';
    i_znak= 0;
    riadok[0]= '\0';
    xyz= x;
    seq= 1;
    line_count=1;

    complet=0;

    prev_char_numberlike= 0;
    curr_char_numberlike= 0;
    numberlike_in_line_counter= 0;
    the_last_coord= 0;
    prev_line_nonblank= 0;
    curr_line_nonblank= 0;


    switch (o_crea_ent) {
    case 'p': {                 /* pts - 116 */
        while (znak != EOF) {
            znak=getc(fr);
            switch (complet) {
            case 0:
                if ((znak >= '0' && znak <= '9') || znak == '-' || znak == '+' ||
                    znak == '.') {
                    sprintf(riadok, "%d,", 116);
                    complet=1;
                    x[i_znak]=znak;
                    i_znak++;
                }
                else {
                    if (znak == ',') {
                        sprintf(riadok, "%d,", 116);
                        complet=1;
                        x[i_znak]='.';
                        i_znak++;
                    }
                }
                break;
            case 1:
                if ((znak >= '0' && znak <= '9') || znak == '-' || znak == '+' ||
                    znak == '.') {
                    x[i_znak]=znak;
                    i_znak++;
                }
                else {
                    switch (znak) {
                    case 'e':
                    case 'E':
                        x[i_znak]='D';
                        i_znak++;
                        break;
                    case ',':
                        x[i_znak]='.';
                        i_znak++;
                        break;
                    case ' ':
                    case '\t':
                        x[i_znak]='\0';
                        if (i_znak != 0) {
                            complet=2;
                            i_znak=0;
                        }
                        break;
                    case '\n':
                        x[i_znak]='\0';
                        strcpy(y, "0D0");
                        strcpy(z, "0D0");
                        sprintf(zbytok, "%s,%s,%s,0;", x, y, z);
                        strcat(riadok, zbytok);
                        pd_linestring_complement(riadok, 2*emark+1, seq);
                        fprintf(fw, "%s\n", riadok);
                        seq++;
                        line_count++;/* condition for complex entities */
                        emark++;
                        complet=0;
                        i_znak=0;
                        break;
                    }
                }
                break;
            case 2:
                if ((znak >= '0' && znak <= '9') || znak == '-' || znak == '+' ||
                    znak == '.') {
                    y[i_znak]=znak;
                    i_znak++;
                }
                else {
                    switch (znak) {
                    case 'e':
                    case 'E':
                        y[i_znak]='D';
                        i_znak++;
                        break;
                    case ',':
                        y[i_znak]='.';
                        i_znak++;
                        break;
                    case ' ':
                    case '\t':
                        y[i_znak]='\0';
                        if (i_znak != 0) {
                            complet=3;
                            i_znak=0;
                        }
                        break;
                    case '\n':
                        y[i_znak]='\0';
                        if (i_znak == 0) strcpy(y, "0D0");
                        strcpy(z, "0D0");
                        sprintf(zbytok, "%s,%s,%s,0;", x, y, z);
                        strcat(riadok, zbytok);
                        pd_linestring_complement(riadok, 2*emark+1, seq);
                        fprintf(fw, "%s\n", riadok);
                        seq++;
                        line_count++;/* condition for complex entities */
                        emark++;
                        complet=0;
                        i_znak=0;
                        break;
                    }
                }
                break;
            case 3:
                if ((znak >= '0' && znak <= '9') || znak == '-' || znak == '+' ||
                    znak == '.') {
                    z[i_znak]=znak;
                    i_znak++;
                }
                else {
                    switch (znak) {
                    case 'e':
                    case 'E':
                        z[i_znak]='D';
                        i_znak++;
                        break;
                    case ',':
                        z[i_znak]='.';
                        i_znak++;
                        break;
                    case ' ':
                    case '\t':
                        z[i_znak]='\0';
                        break;
                    case '\n':
                        z[i_znak]='\0';
                        if (i_znak == 0) strcpy(z, "0D0");
                        sprintf(zbytok, "%s,%s,%s,0;", x, y, z);
                        strcat(riadok, zbytok);
                        pd_linestring_complement(riadok, 2*emark+1, seq);
                        fprintf(fw, "%s\n", riadok);
                        seq++;
                        line_count++;/* condition for complex entities */
                        emark++;
                        complet=0;
                        i_znak=0;
                        break;
                    }
                }
                break;
            }
        }
        break;
    } /*end case p*/
    case 'g': {                 /* groups of pts - 106 */
        while (znak != EOF) {
            znak=getc(fr);
            switch (complet) {
            case 0:
                if ((znak >= '0' && znak <= '9') || znak == '-' || znak == '+' ||
                    znak == '.') {
                    sprintf(riadok, "106,2,%d,", entity[emark][PD_CNT]);
                    complet=1;
                    x[i_znak]=znak;
                    i_znak++;
                }
                else {
                    if (znak == ',') {
                        sprintf(riadok, "106,2,%d,", entity[emark][PD_CNT]);
                        complet=1;
                        x[i_znak]='.';
                        i_znak++;
                    }
                }
                break;
            case 1:
                if ((znak >= '0' && znak <= '9') || znak == '-' || znak == '+' ||
                    znak == '.') {
                    x[i_znak]=znak;
                    i_znak++;
                }
                else {
                    switch (znak) {
                    case 'e':
                    case 'E':
                        x[i_znak]='D';
                        i_znak++;
                        break;
                    case ',':
                        x[i_znak]='.';
                        i_znak++;
                        break;
                    case ' ':
                    case '\t':
                        x[i_znak]='\0';
                        if (i_znak != 0) {
                            complet=2;
                            i_znak=0;
                        }
                        break;
                    case '\n':
                        x[i_znak]='\0';
                        strcpy(y, "0D0");
                        strcpy(z, "0D0");
                        if (line_count == entity[emark][PD_CNT]) {
                            sprintf(zbytok, "%s,%s,%s;", x, y, z);
                            strcat(riadok, zbytok);
                            pd_linestring_complement(riadok, 2*emark+1, seq);
                            fprintf(fw, "%s\n", riadok);
                            complet=0;
                            line_count=1;
                            emark++;
                        }
                        else {
                            sprintf(zbytok, "%s,%s,%s,", x, y, z);
                            strcat(riadok, zbytok);
                            pd_linestring_complement(riadok, 2*emark+1, seq);
                            fprintf(fw, "%s\n", riadok);
                            riadok[0]='\0';
                            line_count++;
                        }
                        seq++;
                        i_znak=0;
                        break;
                    }
                }
                break;
            case 2:
                if ((znak >= '0' && znak <= '9') || znak == '-' || znak == '+' ||
                    znak == '.') {
                    y[i_znak]=znak;
                    i_znak++;
                }
                else {
                    switch (znak) {
                    case 'e':
                    case 'E':
                        y[i_znak]='D';
                        i_znak++;
                        break;
                    case ',':
                        y[i_znak]='.';
                        i_znak++;
                        break;
                    case ' ':
                    case '\t':
                        y[i_znak]='\0';
                        if (i_znak != 0) {
                            complet=3;
                            i_znak=0;
                        }
                        break;
                    case '\n':
                        y[i_znak]='\0';
                        if (i_znak == 0) strcpy(y, "0D0");
                        strcpy(z, "0D0");
                        if (line_count == entity[emark][PD_CNT]) {
                            sprintf(zbytok, "%s,%s,%s;", x, y, z);
                            strcat(riadok, zbytok);
                            pd_linestring_complement(riadok, 2*emark+1, seq);
                            fprintf(fw, "%s\n", riadok);
                            complet=0;
                            line_count=1;
                            emark++;
                        }
                        else {
                            sprintf(zbytok, "%s,%s,%s,", x, y, z);
                            strcat(riadok, zbytok);
                            pd_linestring_complement(riadok, 2*emark+1, seq);
                            fprintf(fw, "%s\n", riadok);
                            complet=1;
                            riadok[0]='\0';
                            line_count++;
                        }
                        seq++;
                        i_znak=0;
                        break;
                    }
                }
                break;
            case 3:
                if ((znak >= '0' && znak <= '9') || znak == '-' || znak == '+' ||
                    znak == '.') {
                    z[i_znak]=znak;
                    i_znak++;
                }
                else {
                    switch (znak) {
                    case 'e':
                    case 'E':
                        z[i_znak]='D';
                        i_znak++;
                        break;
                    case ',':
                        z[i_znak]='.';
                        i_znak++;
                        break;
                    case ' ':
                    case '\t':
                        z[i_znak]='\0';
                        break;
                    case '\n':
                        z[i_znak]='\0';
                        if (i_znak == 0) strcpy(z, "0D0");
                        if (line_count == entity[emark][PD_CNT]) {
                            sprintf(zbytok, "%s,%s,%s;", x, y, z);
                            strcat(riadok, zbytok);
                            pd_linestring_complement(riadok, 2*emark+1, seq);
                            fprintf(fw, "%s\n", riadok);
                            complet=0;
                            line_count=1;
                            emark++;
                        }
                        else {
                            sprintf(zbytok, "%s,%s,%s,", x, y, z);
                            strcat(riadok, zbytok);
                            pd_linestring_complement(riadok, 2*emark+1, seq);
                            fprintf(fw, "%s\n", riadok);
                            riadok[0]='\0';
                            complet=1;
                            line_count++;
                        }
                        seq++;
                        i_znak=0;
                        break;
                    }
                }
                break;
            }
        }
        break;
    } /*end case g*/
    case 'c': {                 /* curves - 126 */
        /* State variables: -global */
        emark= 0;
        fseek(fr, 0L, SEEK_SET);

        /* State variables: -local */
        znak= ' ';
        i_znak= 0;
        riadok[0]= '\0';
        xyz= x;
        seq= 1;
        line_count= 0;                  /* ! - different than for 116 and 106*/
        prev_char_numberlike= 0;
        curr_char_numberlike= 0;
        numberlike_in_line_counter= 0;
        the_last_coord= 0;        /* {NOTHING, X, Y, Z} */
        prev_line_nonblank= 0;
        curr_line_nonblank= 0;


        while (znak != EOF) {
            znak=getc(fr);


            /* At first, proces chars */
            prev_char_numberlike= curr_char_numberlike;
            curr_char_numberlike= is_numberlike(znak);
            if (curr_char_numberlike == 1) {
                numberlike_in_line_counter++;  /* count valid chars */
                /* replace non-standard characters */
                switch (znak) {
                case 'e':
                case 'E':
                case 'd':
                case 'D':
                    znak= 'D';
                    break;
                case ',':
                    znak= '.';
                    break;
                }
            }

            if (prev_char_numberlike == 0 && curr_char_numberlike == 1) {         /* chars pattern 0 1 */
                the_last_coord++;
                i_znak= 0;
                switch (the_last_coord) {
                case 1:
                    xyz=x;
                    break;
                case 2:
                    xyz=y;
                    break;
                case 3:
                    xyz=z;
                    break;
                }
                xyz[i_znak]=znak;
                i_znak++;
            }
            else if (prev_char_numberlike == 1 && curr_char_numberlike == 1) {    /* chars pattern 1 1 */
                xyz[i_znak]=znak;
                i_znak++;
            }
            else if (prev_char_numberlike == 1 && curr_char_numberlike == 0) {    /* chars pattern 1 0 */
                xyz[i_znak]='\0';
                i_znak++;
            }


            /* on \n and EOF */
            if (znak == '\n' || znak == EOF) {
                prev_line_nonblank = curr_line_nonblank;
                if (numberlike_in_line_counter > 0) {           /* Is current line non-blank? */
                    curr_line_nonblank= 1;
                    line_count++;
                }
                else {
                    curr_line_nonblank= 0;
                }
                numberlike_in_line_counter= 0;                    /* already used */
                prev_char_numberlike = 0;                         /* ready for next line */
                curr_char_numberlike = 0;                         /* ready for next line */

                if (prev_line_nonblank == 0 && curr_line_nonblank == 1) {             /* write header */
                    sprintf(riadok, "126,%d,1,0,0,1,0,", entity_coords[emark] -1);                   /* line 1 */
                    pd_linestring_complement(riadok, 2*emark+1, seq);
                    fprintf(fw, "%s\n", riadok);
                    seq++;
                    sprintf(riadok, "0.0D0,0.0D0,");                                              /* knot seq - the first line */
                    pd_linestring_complement(riadok, 2*emark+1, seq);
                    fprintf(fw, "%s\n", riadok);
                    seq++;
                    for (i=1; i < entity_coords[emark] -1; i++) {
                        sprintf(riadok, "%d.0D0,", i);                                              /*      - intermediate lines */
                        pd_linestring_complement(riadok, 2*emark+1, seq);
                        fprintf(fw, "%s\n", riadok);
                        seq++;
                    }
                    sprintf(riadok, "%d.0D0,%d.0D0,", entity_coords[emark] -1, entity_coords[emark] -1);  /*      - the last line */
                    pd_linestring_complement(riadok, 2*emark+1, seq);
                    fprintf(fw, "%s\n", riadok);
                    seq++;
                    for (i=1; i <= entity_coords[emark]; i++) {
                        sprintf(riadok, "%d.0D0,", i);   /* weights */
                        pd_linestring_complement(riadok, 2*emark+1, seq);
                        fprintf(fw, "%s\n", riadok);
                        seq++;
                    }
                }

                if (curr_line_nonblank == 1) {                                        /* write xyz */
                    switch (the_last_coord) {
                    case 1:
                        strcpy(y, "0D0");
                        strcpy(z, "0D0");
                        break;
                    case 2:
                        strcpy(z, "0D0");
                        break;
                    case 3:
                        break;
                    }
                    sprintf(riadok, "%s,%s,%s,", x, y, z);
                    pd_linestring_complement(riadok, 2*emark+1, seq);
                    fprintf(fw, "%s\n", riadok);
                    seq++;
                    the_last_coord= 0;
                    i_znak= 0;
                }

                if ((prev_line_nonblank == 1 && curr_line_nonblank == 0)
                    || (znak == EOF && curr_line_nonblank == 1) ) {                   /* write tail */
                    sprintf(riadok, "0.0D0,%d.0D0,0.0D0,0.0D0,0.0D0;", entity_coords[emark] -1 -1);                   /* line 1 */
                    pd_linestring_complement(riadok, 2*emark+1, seq);
                    fprintf(fw, "%s\n", riadok);
                    seq++;
                    emark++;
                }
            } /* end if (znak == '\n' || znak == EOF) */
        } /* end while (znak != EOF) */
        break;
    } /* end case c */
    case 's': {                 /* surface - 128 */
        /* State variables: -global */
        emark= 0;
        fseek(fr, 0L, SEEK_SET);

        /* State variables: -local */
        znak= ' ';
        i_znak= 0;
        riadok[0]= '\0';
        xyz= x;
        seq= 1;
        line_count= 0;                  /* ! - different than for 116 and 106*/
        prev_char_numberlike= 0;
        curr_char_numberlike= 0;
        numberlike_in_line_counter= 0;
        the_last_coord= 0;        /* {NOTHING, X, Y, Z} */
        prev_line_nonblank= 0;
        curr_line_nonblank= 0;


        while (znak != EOF) {
            znak=getc(fr);


            /* At first, proces chars */
            prev_char_numberlike= curr_char_numberlike;
            curr_char_numberlike= is_numberlike(znak);
            if (curr_char_numberlike == 1) {
                numberlike_in_line_counter++;  /* count valid chars */
                /* replace non-standard characters */
                switch (znak) {
                case 'e':
                case 'E':
                case 'd':
                case 'D':
                    znak= 'D';
                    break;
                case ',':
                    znak= '.';
                    break;
                }
            }

            if (prev_char_numberlike == 0 && curr_char_numberlike == 1) {         /* chars pattern 0 1 */
                the_last_coord++;
                i_znak= 0;
                switch (the_last_coord) {
                case 1:
                    xyz=x;
                    break;
                case 2:
                    xyz=y;
                    break;
                case 3:
                    xyz=z;
                    break;
                }
                xyz[i_znak]=znak;
                i_znak++;
            }
            else if (prev_char_numberlike == 1 && curr_char_numberlike == 1) {    /* chars pattern 1 1 */
                xyz[i_znak]=znak;
                i_znak++;
            }
            else if (prev_char_numberlike == 1 && curr_char_numberlike == 0) {    /* chars pattern 1 0 */
                xyz[i_znak]='\0';
                i_znak++;
            }


            /* on \n and EOF */
            if (znak == '\n' || znak == EOF) {
                prev_line_nonblank = curr_line_nonblank;
                if (numberlike_in_line_counter > 0) {           /* Is current line non-blank? */
                    curr_line_nonblank= 1;
                    line_count++;
                }
                else {
                    curr_line_nonblank= 0;
                }
                numberlike_in_line_counter= 0;                    /* already used */
                prev_char_numberlike = 0;                         /* ready for next line */
                curr_char_numberlike = 0;                         /* ready for next line */

                if (prev_line_nonblank == 0 && curr_line_nonblank == 1) {             /* write header */
                    sprintf(riadok, "126,%d,1,0,0,1,0,", entity_coords[emark] -1);                   /* line 1 */
                    pd_linestring_complement(riadok, 2*emark+1, seq);
                    fprintf(fw, "%s\n", riadok);
                    seq++;
                    sprintf(riadok, "0.0D0,0.0D0,");                                              /* knot seq - the first line */
                    pd_linestring_complement(riadok, 2*emark+1, seq);
                    fprintf(fw, "%s\n", riadok);
                    seq++;
                    for (i=1; i < entity_coords[emark] -1; i++) {
                        sprintf(riadok, "%d.0D0,", i);                                              /*      - intermediate lines */
                        pd_linestring_complement(riadok, 2*emark+1, seq);
                        fprintf(fw, "%s\n", riadok);
                        seq++;
                    }
                    sprintf(riadok, "%d.0D0,%d.0D0,", entity_coords[emark] -1, entity_coords[emark] -1);  /*      - the last line */
                    pd_linestring_complement(riadok, 2*emark+1, seq);
                    fprintf(fw, "%s\n", riadok);
                    seq++;
                    for (i=1; i <= entity_coords[emark]; i++) {
                        sprintf(riadok, "%d.0D0,", i);   /* weights */
                        pd_linestring_complement(riadok, 2*emark+1, seq);
                        fprintf(fw, "%s\n", riadok);
                        seq++;
                    }
                }

                if (curr_line_nonblank == 1) {                                        /* write xyz */
                    switch (the_last_coord) {
                    case 1:
                        strcpy(y, "0D0");
                        strcpy(z, "0D0");
                        break;
                    case 2:
                        strcpy(z, "0D0");
                        break;
                    case 3:
                        break;
                    }
                    sprintf(riadok, "%s,%s,%s,", x, y, z);
                    pd_linestring_complement(riadok, 2*emark+1, seq);
                    fprintf(fw, "%s\n", riadok);
                    seq++;
                    the_last_coord= 0;
                    i_znak= 0;
                }

                if ((prev_line_nonblank == 1 && curr_line_nonblank == 0)
                    || (znak == EOF && curr_line_nonblank == 1) ) {                   /* write tail */
                    sprintf(riadok, "0.0D0,%d.0D0,0.0D0,0.0D0,0.0D0;", entity_coords[emark] -1 -1);                   /* line 1 */
                    pd_linestring_complement(riadok, 2*emark+1, seq);
                    fprintf(fw, "%s\n", riadok);
                    seq++;
                    emark++;
                }
            } /* end if (znak == '\n' || znak == EOF) */
        } /* end while (znak != EOF) */
        break;
    } /* end case s */

    } /* end switch (o_crea_ent) */
}

void writeiges_t(int ** entity, int & entity_sum, FILE * fw)
{
    /* write  T section to IGES file */
    char riadok[81];
    char retazec[9];
    int pd;
    /*1      81      81      81      81      81      81      81      81      81      8*/
    strcpy(riadok,   "S0000001G0000006D000000xP000000x                                        T      1");
    sprintf(retazec, "%07d", (2*entity_sum));
    fill_s(retazec, 17, 7, riadok);
    pd=entity[entity_sum-1][PD_PTR]+entity[entity_sum-1][PD_CNT]-1;
    sprintf(retazec, "%07d", pd);
    fill_s(retazec, 25, 7, riadok);
    fprintf(fw, "%s\n", riadok);
}

int build_de_c(int ** entity, int & emark, int & entity_sum)
{
    /* Build IGES DE section in memory - create curve - 126*/
    int  current_line_chars,    /* to detect non-empty line (number of non-blank characters in line) */
        previous_line_chars,   /* to detect that the current line is not the first (number of non-blank characters in previous line) */
        coordinates,           /* number of coordinates (contiguous lines) */
        param_seq;             /* pointer to PD data */
    char character;
    static int lay;
    current_line_chars=0;
    previous_line_chars=0;
    coordinates=0;
    param_seq=1;
    lay=0;
    character=' ';


    FILE * fr =NULL;//REMOVE (input txt file)
    fseek(fr, 0L, SEEK_SET);
    do {
        character=getc(fr);
        switch (character) {
        case '\n':
            if (current_line_chars != 0) {                       /* \n after non-empty line (with coordinates) */
                coordinates++;
            }
            else if (previous_line_chars != 0) {                 /* \n in empty line after contiguous block */
                entity_coords[emark]=coordinates;
                entity[emark][E_TYPE]=126;
                entity[emark][PD_PTR]=param_seq;
                entity[emark][PD_CNT]=2+3*coordinates;
                entity[emark][ELAYER]=lay;
                param_seq+=entity[emark][PD_CNT];
                coordinates=0;
                entity_sum++;
                if (entity_sum == MAXENTITY) return -1;
                emark++;
                if (o_incr_layer == 1) lay++;
            }
            previous_line_chars=current_line_chars;
            current_line_chars=0;
            break;
        case EOF:
            if (previous_line_chars != 0) {
                entity_coords[emark]=coordinates;
                entity[emark][E_TYPE]=126;
                entity[emark][PD_PTR]=param_seq;
                entity[emark][PD_CNT]=2+3*coordinates;
                entity[emark][ELAYER]=lay;
                param_seq+=entity[emark][PD_CNT];
                coordinates=0;
                entity_sum++;
                if (entity_sum == MAXENTITY) return -1;
                emark++;
                if (o_incr_layer == 1) lay++;
            }
            break;
        case ' ':
        case '\t':
        case '\r': break;
        default: current_line_chars++; break;            /* detect non-empty line */
        }
    } while (character != EOF);
    return 0;
}

int build_de_s(int ** entity, int & emark, int & entity_sum)
{
    /* Build IGES DE section in memory - create surface - 128*/
    int  current_line_chars,    /* to detect non-empty line (number of non-blank characters in line) */
        previous_line_chars,   /* to detect that the current line is not the first (number of non-blank characters in previous line) */
        coordinates,           /* number of coordinates (contiguous lines) */
        param_seq;             /* pointer to PD data */
    char character;
    static int lay;
    current_line_chars=0;
    previous_line_chars=0;
    coordinates=0;
    param_seq=1;
    lay=0;
    character=' ';
    FILE * fr =NULL;//REMOVE (input txt file)
    fseek(fr, 0L, SEEK_SET);
    do {
        character=getc(fr);
        switch (character) {
        case '\n':
            if (current_line_chars != 0) {  /* \n after non-empty line (with coordinates) */
                coordinates++;
            }
            else if (previous_line_chars != 0) {  /* \n in empty line after contiguous block */
                entity_coords[emark]=coordinates;
                entity[emark][E_TYPE]=128;
                entity[emark][PD_PTR]=param_seq;
                entity[emark][PD_CNT]=2+3*coordinates;
                entity[emark][ELAYER]=lay;
                param_seq+=entity[emark][PD_CNT];
                coordinates=0;
                entity_sum++;
                if (entity_sum == MAXENTITY) return -1;
                emark++;
                if (o_incr_layer == 1) lay++;
            }
            previous_line_chars=current_line_chars;
            current_line_chars=0;
            break;
        case EOF:
            if (previous_line_chars != 0) {
                entity_coords[emark]=coordinates;
                entity[emark][E_TYPE]=128;
                entity[emark][PD_PTR]=param_seq;
                entity[emark][PD_CNT]=2+3*coordinates;
                entity[emark][ELAYER]=lay;
                param_seq+=entity[emark][PD_CNT];
                coordinates=0;
                entity_sum++;
                if (entity_sum == MAXENTITY) return -1;
                emark++;
                if (o_incr_layer == 1) lay++;
            }
            break;
        case ' ':
        case '\t':
        case '\r': break;
        default: current_line_chars++; break;  /* detect non-empty line */
        }
    } while (character != EOF);
    return 0;
}



template<class T>
void gsFileData<T>::writeIges(String const & fname)
{    
    FILE * fw = fopen((fname+".igs").c_str(), "w");
    int (*build_de)(void); /* pointer to function writing the DE */

    /*   1      81      81      81      81      81      81      81      8*/
    fprintf(fw,"IGES file Created by G+Smo (https://github.com/gismo/gismo/wiki)"
            "        S      1\n");
    fprintf(fw,"1H,,1H;,,7Hout.txt,                                             "
            "        G      1\n");
    fprintf(fw,"19Hpegas freeware 2010,12HIGSPREAD 0.5,32,38,6,308,15,,         "
            "        G      2\n");
    fprintf(fw,"1.0D0,2,2HMM,1000,0.1D0,15H20001112.123000,                     "
            "        G      3\n");
    fprintf(fw,"0.000000000000001D0,16D0,                                       "
            "        G      4\n");
    fprintf(fw,"16HPeter Gasparovic,17HAir Force Academy,11,0,                  "
            "        G      5\n");
    fprintf(fw,"15H20001112.123000;                                             "
            "        G      6\n");
    switch (o_crea_ent)
    {
//  case 'p': build_de= build_de_p; break;
//  case 'g': build_de= build_de_g; break;
    case 'c': build_de= build_de_c; break;
    case 's': build_de= build_de_s; break;
    }
    int entity_sum=0;
    int emark=0;
    GISMO_ENSURE(build_de() == -1,//run build_de() function
                 "Insufficient limit \"MAXENTITY\" for entity processing.");
    writeiges_de(fw);
    writeiges_pd(fw);
    writeiges_t(fw);

    const bool fc = (fclose(fw) != EOF);
    if (fc) gsWarn<< "Error: File closing didn't succeeded!\n";
}

