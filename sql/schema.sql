create table if not exists images (
    id integer primary key autoincrement,
    file_path text not null,
    filter text not null,
    time text not null,
    asteroid_value real not null,
    calibrated_value real not null,
    ref_star_name text not null,

    foreign key (ref_star_name) references ref_stars(name)
);

create table if not exists ref_stars (
    id integer primary key autoincrement,
    name text not null unique,
    abs_mag real not null,
    distance real not null,
    img_num integer not null
);

