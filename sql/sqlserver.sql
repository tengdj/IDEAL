create table polygon (id int not null, geog geography, constraint PK_T PRIMARY KEY (ID))
create table polygon_small (id int not null, geog geography, constraint PK_T_small PRIMARY KEY (ID))

create table point (id int not null, geog geography, constraint PK_T PRIMARY KEY (ID))
CREATE SPATIAL INDEX polygon_grid on polygon(geog) USING GEOGRAPHY_AUTO_GRID
select count(*) from point t1, polygon t2 where t2.geog.STContains(t1.geog)=1 and t1.id=1

CREATE SPATIAL INDEX polygon_grid ON polygon(geog) USING GEOGRAPHY_GRID WITH ( BOUNDING_BOX = ( xmin=-180, ymin=-90, xmax=180, ymax=90), GRIDS = (MEDIUM, LOW, MEDIUM, HIGH ), CELLS_PER_OBJECT = 640, PAD_INDEX  = ON );  

CREATE SPATIAL INDEX grid_index on polygon_small(geog) USING GEOGRAPHY_AUTO_GRID

CREATE SPATIAL INDEX grid_index ON polygon_small(geog) USING GEOGRAPHY_GRID WITH ( GRIDS = (HIGH, LOW, low, HIGH ), CELLS_PER_OBJECT = 640, PAD_INDEX  = ON ); 

CREATE SPATIAL INDEX polygon_grid ON polygon(geog) USING GEOGRAPHY_GRID WITH ( BOUNDING_BOX = ( xmin=-180, ymin=-90, xmax=180, ymax=90), GRIDS = (MEDIUM, LOW, MEDIUM, HIGH ), CELLS_PER_OBJECT = 640, PAD_INDEX  = ON );  


time sqlcmd -S localhost -U SA -P '$Terry08161043' -d tengdb -I -Q "select count(*) from point t1, polygon_small t2 where t2.geog.STContains(t1.geog)=1;"
sqlcmd -S localhost -U SA -P '$Terry08161043' -d tengdb -I -i polygon.small.sql