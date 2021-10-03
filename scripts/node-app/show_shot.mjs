
import axios from 'axios';
import fs from 'fs';


const graphicsClient = axios.create({baseURL: 'http://127.0.0.1:18082/', timeout: 1000});
const projectorClient = axios.create({baseURL: 'http://127.0.0.1:18080/', timeout: 1000});


const getGraphics = async (table, locations, shotInfo) => {
    const {data: {message, graphics, success}, status} = await graphicsClient.post(
		'/shot-info/',
        {
            params: {
                table,
                locations,
                'shot-info': shotInfo,
            },
        });
    if (!success || status !== 200) {
        console.log('Unable to get graphics', message);
    }
    return graphics;
};


const pushGraphics = async graphics => {
    const {data: {message, success}, status} = await projectorClient.put(
		'/graphics/',
		{graphics}
	);
    if (!success || status !== 200) {
        console.log('Unable to push graphics', message);
    }
};

const main = async () => {
	const table = JSON.parse(fs.readFileSync('inputs/table.json'));
	const locations = JSON.parse(fs.readFileSync('inputs/locations.json'));
	const shotInfo = JSON.parse(fs.readFileSync('inputs/shot_info.json'));
	
	console.log('shotInfo', JSON.stringify(shotInfo, null, 2));
	
	const graphics = await getGraphics(table, locations, shotInfo);
	await pushGraphics(graphics);
}

main();

