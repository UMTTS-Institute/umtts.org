<?php
declare(strict_types=1);

header('Content-Type: application/json; charset=utf-8');
header('Cache-Control: no-store, no-cache, must-revalidate, max-age=0');

$cookieName = 'umtts_visitor_id';
$oneYear = time() + 31536000;

if (empty($_COOKIE[$cookieName]) || !preg_match('/^[a-f0-9]{32}$/', $_COOKIE[$cookieName])) {
    $visitorId = bin2hex(random_bytes(16));
    setcookie($cookieName, $visitorId, [
        'expires' => $oneYear,
        'path' => '/',
        'secure' => (!empty($_SERVER['HTTPS']) && $_SERVER['HTTPS'] !== 'off'),
        'httponly' => true,
        'samesite' => 'Lax'
    ]);
} else {
    $visitorId = $_COOKIE[$cookieName];
}

$page = $_GET['page'] ?? '/';
$page = preg_replace('/[#?].*$/', '', (string)$page);
$page = preg_replace('/[^a-zA-Z0-9_\/\.\-]/', '', $page);
$page = $page === '' ? '/' : $page;

$dir = __DIR__ . DIRECTORY_SEPARATOR . 'data';
$file = $dir . DIRECTORY_SEPARATOR . 'visitor-counts.json';

if (!is_dir($dir)) {
    mkdir($dir, 0755, true);
}

$visitorHash = hash('sha256', $visitorId);
$now = gmdate('c');

$fp = fopen($file, 'c+');
if (!$fp) {
    http_response_code(500);
    echo json_encode(['error' => 'counter storage unavailable']);
    exit;
}

flock($fp, LOCK_EX);
$raw = stream_get_contents($fp);
$data = json_decode($raw ?: '', true);

if (!is_array($data)) {
    $data = [
        'site' => ['visits' => 0, 'visitors' => []],
        'pages' => [],
        'updated_utc' => $now
    ];
}

if (!isset($data['site']['visits'])) {
    $data['site']['visits'] = 0;
}
if (!isset($data['site']['visitors']) || !is_array($data['site']['visitors'])) {
    $data['site']['visitors'] = [];
}
if (!isset($data['pages'][$page])) {
    $data['pages'][$page] = ['visits' => 0, 'visitors' => []];
}
if (!isset($data['pages'][$page]['visitors']) || !is_array($data['pages'][$page]['visitors'])) {
    $data['pages'][$page]['visitors'] = [];
}

$data['site']['visits']++;
$data['site']['visitors'][$visitorHash] = $now;

$data['pages'][$page]['visits']++;
$data['pages'][$page]['visitors'][$visitorHash] = $now;

$data['updated_utc'] = $now;

rewind($fp);
ftruncate($fp, 0);
fwrite($fp, json_encode($data, JSON_PRETTY_PRINT | JSON_UNESCAPED_SLASHES));
fflush($fp);
flock($fp, LOCK_UN);
fclose($fp);

echo json_encode([
    'page' => $page,
    'site_visitors' => count($data['site']['visitors']),
    'page_visitors' => count($data['pages'][$page]['visitors']),
    'site_visits' => $data['site']['visits'],
    'page_visits' => $data['pages'][$page]['visits'],
    'updated_utc' => $now
]);
?>
