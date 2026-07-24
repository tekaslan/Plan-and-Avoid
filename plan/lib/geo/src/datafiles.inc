/* 
 * @file
 * @brief Files to read Lat/Lon data from Files
*/







bool geo_readgga(const char * const line, struct Pos * const pos)
{
    const char * aux = line;
    
    /* Check if it is a correct line */
    if(strncmp(line,"$GPGGA",5))   return 1;

    /* Check if it is a correct fix */
    for(uint_fast8_t i = 0u; i < 6u; i++){
        aux = strchr(++aux,',');
        if(!aux) return 1;
    }
    
    /* Get data */
    aux = strchr(line,',');
    aux = strchr(++aux,',');

    pos->lat = floor(atof(++aux) / 100.0);
    pos->lat += (atof(aux) - 100.0 * pos->lat) / 60.0;

    aux = strchr(aux,',');
    if(aux[1] == 'S') pos->lat *= -1.0;

    aux = strchr(++aux,',');
    
    pos->lon = floor(atof(++aux) / 100.0);                                                                                     
    pos->lon += (atof(aux) - 100.0 * pos->lon) / 60.0; 
    
    aux = strchr(aux,',');
    if(aux[1] == 'W') pos->lon *= -1.0;

    /* Read altitude */
    for(uint_fast8_t i = 0u; i < 4u; i++){
        aux = strchr(++aux,',');
        if(!aux) return 1;
    }
    pos->alt = atof(++aux);
    
    /* Checks altitude in meters */
    aux = strchr(aux,',');
    if(++aux != 'M') return 1;
    else             pos->alt *= M_2_FT;

    return 0;
}

bool geo_readcsv(const char * const line, struct Pos * const pos,
                 const uint_fast8_t const cols[3], const char delimiter )
{

    for(uint_fast8_t i = 0; i < 3; i++){
        const char * aux = line;
        if(col[i] > 0){
            for(uint_fast8_t j = 0; j < cols[i]; j++){
                aux = strchr(++aux, delimiter);
            }
            switch(i){
            case 0:     pos->lat = atof(++aux);
            case 1:     pos->lon = atof(++aux);
            case 2:     pos->alt = atof(++aux);
            }
        }else{
            switch(i){
            case 0:     pos->lon = 91.0;
            case 1:     pos->lat = 181.0;
            case 2:     pos->alt = -1.0e5;
            }
        }
    }
    
    return 0;
}

